"""biofile MCP server.

Exposes the functions in the biofile package as MCP tools. Claude reads each
function's name, type hints and docstring to decide when to call it and with
which arguments, so this file contains no text-to-command logic.

Run:  python server.py

Settings (environment variables):
    BIOFILE_ALLOWED_DIR  only files inside this folder can be read (default: current folder)
    BIOFILE_LOG_LEVEL    DEBUG, INFO (default), WARNING or ERROR
    BIOFILE_LOG_FILE     also write logs to this file

Logs go to stderr, never stdout: over stdio, stdout carries the MCP messages
themselves, and any extra output there would break the connection.
"""
import logging
import os
import sys
import time
from functools import wraps
from pathlib import Path

from mcp.server.mcpserver import MCPServer  # MCP Python SDK v2 (was FastMCP in v1)
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from biofile import bam_stats, fastq_stats, h5ad_summary, vcf_stats

ALLOWED_DIR = Path(os.environ.get("BIOFILE_ALLOWED_DIR", Path.cwd())).expanduser().resolve()

# ---- Logging ---------------------------------------------------------------
handlers = [logging.StreamHandler(sys.stderr)]
if os.environ.get("BIOFILE_LOG_FILE"):
    handlers.append(logging.FileHandler(os.environ["BIOFILE_LOG_FILE"]))

logging.basicConfig(
    level=os.environ.get("BIOFILE_LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=handlers,
)
log = logging.getLogger("biofile")


def resolve_path(path: str) -> Path:
    """Turn a user/model-supplied path into an absolute path inside ALLOWED_DIR."""
    full = Path(path).expanduser()
    if not full.is_absolute():
        full = ALLOWED_DIR / full
    full = full.resolve()

    if full != ALLOWED_DIR and ALLOWED_DIR not in full.parents:
        raise ToolError(f"{full} is outside the allowed folder {ALLOWED_DIR}")
    if not full.exists():
        raise ToolError(f"{full} does not exist")
    return full


def only_allowed_folder(fn):
    """Check the `path` argument, log the call, and pass errors back to Claude."""

    @wraps(fn)  # keeps the original name, docstring and arguments for Claude
    def wrapper(path: str, *args, **kwargs):
        start = time.perf_counter()
        log.info("call %s path=%s args=%s", fn.__name__, path, kwargs)

        try:
            full = resolve_path(path)
            result = fn(str(full), *args, **kwargs)
        except ToolError as e:
            log.warning("refused %s: %s", fn.__name__, e)
            raise
        except Exception as e:
            log.exception("failed %s after %.2fs", fn.__name__, time.perf_counter() - start)
            # MCP v2 hides normal Python errors from Claude. ToolError passes the
            # message through, so Claude can see what went wrong and fix its call.
            raise ToolError(f"{type(e).__name__}: {e}") from e

        log.info("done %s in %.2fs", fn.__name__, time.perf_counter() - start)
        return result

    return wrapper


mcp = MCPServer("biofile")

# readOnlyHint tells the client these tools never change anything
read_only = ToolAnnotations(readOnlyHint=True)

TOOLS = [vcf_stats, bam_stats, h5ad_summary, fastq_stats]
for tool in TOOLS:
    mcp.tool(annotations=read_only)(only_allowed_folder(tool))


if __name__ == "__main__":
    log.info("starting biofile server; allowed folder: %s", ALLOWED_DIR)
    mcp.run()  # talks to Claude over stdin/stdout
