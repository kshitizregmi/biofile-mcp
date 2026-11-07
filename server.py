"""biofile MCP server.

Exposes the functions in the biofile package as MCP tools. Claude reads each
function's name, type hints and docstring to decide when to call it and with
which arguments, so this file contains no text-to-command logic.

Run:  python server.py
Only files inside BIOFILE_ALLOWED_DIR (default: the current folder) can be read.
"""
import os
from functools import wraps
from pathlib import Path

from mcp.server.mcpserver import MCPServer  # MCP Python SDK v2 (was FastMCP in v1)
from mcp.server.mcpserver.exceptions import ToolError
from mcp.types import ToolAnnotations

from biofile import bam_stats, fastq_stats, h5ad_summary, vcf_stats

ALLOWED_DIR = Path(os.environ.get("BIOFILE_ALLOWED_DIR", Path.cwd())).expanduser().resolve()


def only_allowed_folder(fn):
    """Check the `path` argument before calling the real function."""

    @wraps(fn)  # keeps the original name, docstring and arguments for Claude
    def wrapper(path: str, *args, **kwargs):
        full = Path(path).expanduser()
        if not full.is_absolute():
            full = ALLOWED_DIR / full
        full = full.resolve()

        if full != ALLOWED_DIR and ALLOWED_DIR not in full.parents:
            raise ToolError(f"{full} is outside the allowed folder {ALLOWED_DIR}")
        if not full.exists():
            raise ToolError(f"{full} does not exist")

        try:
            return fn(str(full), *args, **kwargs)
        except Exception as e:
            # MCP v2 hides normal Python errors from Claude. ToolError passes the
            # message through, so Claude can see what went wrong and fix its call.
            raise ToolError(f"{type(e).__name__}: {e}") from e

    return wrapper


mcp = MCPServer("biofile")

# readOnlyHint tells the client these tools never change anything
read_only = ToolAnnotations(readOnlyHint=True)

for tool in [vcf_stats, bam_stats, h5ad_summary, fastq_stats]:
    mcp.tool(annotations=read_only)(only_allowed_folder(tool))


if __name__ == "__main__":
    mcp.run()  # talks to Claude over stdin/stdout
