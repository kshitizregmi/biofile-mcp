#!/usr/bin/env bash
# Set up biofile-mcp in one step:
#   1. create (or update) the "biofile" conda environment
#   2. make test files and run the test suite
#   3. connect the server to Claude Code
#
# Usage:
#   bash setup.sh                      # server may read ./test_data
#   bash setup.sh /path/to/your/data   # server may read that folder

set -euo pipefail

ENV_NAME="biofile"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${1:-$PROJECT_DIR/test_data}"

cd "$PROJECT_DIR"

# 1. Find conda
if ! command -v conda >/dev/null 2>&1; then
    for base in "$HOME/miniconda3" "$HOME/miniconda" "$HOME/anaconda3" "$HOME/miniforge3"; do
        if [[ -f "$base/etc/profile.d/conda.sh" ]]; then
            source "$base/etc/profile.d/conda.sh"
            break
        fi
    done
fi
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found. Install miniconda first."; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

# 2. Create or update the environment
if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo ">> Updating conda env '$ENV_NAME'"
    conda env update -n "$ENV_NAME" -f environment.yml --prune
else
    echo ">> Creating conda env '$ENV_NAME'"
    conda env create -n "$ENV_NAME" -f environment.yml
fi

conda activate "$ENV_NAME"
ENV_PYTHON="$(command -v python)"
echo ">> Using Python: $ENV_PYTHON"

# 3. Run the tests
echo ">> Making test files"
python make_test_data.py

echo ">> Running the test suite"
python -m pytest -q
echo "   all tests passed"

# 4. Register the server with Claude Code
mkdir -p "$DATA_DIR"
DATA_DIR="$(cd "$DATA_DIR" && pwd)"

if command -v claude >/dev/null 2>&1; then
    echo ">> Registering with Claude Code (allowed folder: $DATA_DIR)"
    claude mcp remove biofile --scope user >/dev/null 2>&1 || true
    claude mcp add biofile --scope user \
        -e BIOFILE_ALLOWED_DIR="$DATA_DIR" \
        -- "$ENV_PYTHON" "$PROJECT_DIR/server.py"
    claude mcp list
else
    echo ">> 'claude' command not found. Register it yourself later with:"
    echo "   claude mcp add biofile --scope user -e BIOFILE_ALLOWED_DIR=\"$DATA_DIR\" -- \"$ENV_PYTHON\" \"$PROJECT_DIR/server.py\""
fi

echo
echo "Done. Open Claude Code and ask, for example:"
echo "   \"What's the Ts/Tv ratio of test.vcf?\""
