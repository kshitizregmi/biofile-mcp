#!/bin/sh
# Builds the Python environment the biofile server runs in.
#
# The environment lives in the plugin data folder, so it survives plugin
# updates. It is rebuilt only when requirements-server.txt changes, which
# makes every run after the first one finish straight away.
#
# All output goes to stderr. Claude Code adds a session hook's stdout to the
# conversation, and the MCP server uses stdout for its protocol.
set -e

ROOT="${CLAUDE_PLUGIN_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
DATA="${CLAUDE_PLUGIN_DATA:-$HOME/.biofile-mcp}"
VENV="$DATA/venv"
REQUIREMENTS="$ROOT/requirements-server.txt"
LOCK="$DATA/install.lock"

mkdir -p "$DATA"

# The session hook and the server can both start on the first run.
# Only one of them installs; the other waits here until it is done.
waited=0
until mkdir "$LOCK" 2>/dev/null; do
    sleep 2
    waited=$((waited + 2))
    if [ "$waited" -ge 600 ]; then
        echo "biofile: clearing an install lock left by an earlier run" >&2
        rmdir "$LOCK" 2>/dev/null || true
        waited=0
    fi
done
trap 'rmdir "$LOCK"' EXIT

# Same requirements as last time: nothing to do.
if cmp -s "$REQUIREMENTS" "$VENV/installed-requirements.txt"; then
    exit 0
fi

# The MCP SDK needs Python 3.10 or newer. Use the newest one we can find.
PYTHON=""
for candidate in python3.13 python3.12 python3.11 python3.10 python3; do
    if command -v "$candidate" >/dev/null 2>&1 &&
        "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 10))' 2>/dev/null; then
        PYTHON="$candidate"
        break
    fi
done
if [ -z "$PYTHON" ]; then
    echo "biofile: Python 3.10 or newer is required." >&2
    echo "biofile: install it (for example 'brew install python@3.12') and restart Claude Code." >&2
    exit 1
fi

echo "biofile: installing Python packages into $VENV (first run only)" >&2
rm -rf "$VENV"
"$PYTHON" -m venv "$VENV" >&2
"$VENV/bin/python" -m pip install --quiet --upgrade pip >&2
"$VENV/bin/python" -m pip install --quiet -r "$REQUIREMENTS" >&2
cp "$REQUIREMENTS" "$VENV/installed-requirements.txt"
echo "biofile: Python packages installed" >&2
