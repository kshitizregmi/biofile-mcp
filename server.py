"""MCP server exposing the biofile summaries as read-only tools.

Each function in the ``biofile`` package is registered as a tool. MCP clients
discover the tools at startup and receive each function's name, parameters and
docstring as the tool's schema and description.

Usage::

    python server.py

Environment variables:

``BIOFILE_ALLOWED_DIR``
    Folder the tools may read from. Relative paths are resolved against it and
    paths outside it are refused. Defaults to the current working directory.
``BIOFILE_LOG_LEVEL``
    ``DEBUG``, ``INFO`` (default), ``WARNING`` or ``ERROR``.
``BIOFILE_LOG_FILE``
    Optional path of a file to write logs to, in addition to stderr.
"""
import logging
import os
import sys
import time
from functools import wraps
from pathlib import Path

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from biofile import bam_stats, fastq_stats, h5ad_summary, vcf_stats

ALLOWED_DIR = Path(os.environ.get("BIOFILE_ALLOWED_DIR", Path.cwd())).expanduser().resolve()

# Logs go to stderr only. With the stdio transport, stdout carries the MCP
# messages, and anything else printed there would corrupt the connection.
_handlers = [logging.StreamHandler(sys.stderr)]
if os.environ.get("BIOFILE_LOG_FILE"):
    _handlers.append(logging.FileHandler(os.environ["BIOFILE_LOG_FILE"]))

logging.basicConfig(
    level=os.environ.get("BIOFILE_LOG_LEVEL", "INFO").upper(),
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=_handlers,
)
log = logging.getLogger("biofile")


def resolve_path(path: str) -> Path:
    """Resolve a requested path and confirm it is a readable location.

    Parameters
    ----------
    path : str
        Absolute path, or a path relative to ``ALLOWED_DIR``. ``~`` is expanded.

    Returns
    -------
    pathlib.Path
        The absolute, symlink-resolved path.

    Raises
    ------
    ToolError
        If the resolved path is outside ``ALLOWED_DIR`` (including via ``..``
        or symlinks), or if it does not exist.
    """
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
    """Wrap a biofile function for use as an MCP tool.

    The wrapper resolves and checks the ``path`` argument with
    :func:`resolve_path`, logs each call with its duration and outcome, and
    converts any exception into a :class:`ToolError`. MCP servers hide the text
    of other exceptions from the client, so this conversion is what lets the
    caller see the actual reason for a failure and correct its request.

    The wrapped function keeps the original name, signature and docstring,
    which become the tool's name, input schema and description.
    """

    @wraps(fn)
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
            raise ToolError(f"{type(e).__name__}: {e}") from e

        log.info("done %s in %.2fs", fn.__name__, time.perf_counter() - start)
        return result

    return wrapper


mcp = MCPServer("biofile")

# Mark every tool read-only so clients know a call never changes anything
read_only = ToolAnnotations(readOnlyHint=True)

TOOLS = [vcf_stats, bam_stats, h5ad_summary, fastq_stats]
for tool in TOOLS:
    mcp.tool(annotations=read_only)(only_allowed_folder(tool))


if __name__ == "__main__":
    log.info("starting biofile server; allowed folder: %s", ALLOWED_DIR)
    mcp.run()
