#!/bin/sh
# Starts the biofile MCP server inside the plugin's Python environment.
#
# Claude Code runs this when the plugin is enabled. If the environment is
# not ready yet, install.sh builds it (or waits for the session hook that
# is already building it) before the server starts.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
DATA="${CLAUDE_PLUGIN_DATA:-$HOME/.biofile-mcp}"

# stdout belongs to the MCP protocol, so install messages go to stderr.
sh "$HERE/install.sh" >&2

exec "$DATA/venv/bin/python" "$HERE/../server.py"
