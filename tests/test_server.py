"""Tests for the MCP server layer.

Covers path checks, error reporting and logging in the tool wrapper, plus an
end-to-end test that starts ``server.py`` as a subprocess and talks to it over
stdio with the MCP client library.
"""
import asyncio
import json
import logging
import os
import sys
from pathlib import Path

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.server.mcpserver.exceptions import ToolError

import server

PROJECT_DIR = Path(__file__).resolve().parents[1]


@pytest.fixture
def allowed(data_dir, monkeypatch):
    monkeypatch.setattr(server, "ALLOWED_DIR", data_dir.resolve())
    return data_dir.resolve()


# Path checks

def test_relative_path_resolves_inside_allowed_folder(allowed):
    assert server.resolve_path("test.vcf") == allowed / "test.vcf"


def test_path_outside_allowed_folder_is_refused(allowed):
    with pytest.raises(ToolError, match="outside the allowed folder"):
        server.resolve_path("/etc/passwd")


def test_dot_dot_cannot_escape(allowed):
    with pytest.raises(ToolError, match="outside the allowed folder"):
        server.resolve_path("../../etc/passwd")


def test_missing_file_is_reported(allowed):
    with pytest.raises(ToolError, match="does not exist"):
        server.resolve_path("nope.vcf")


# Error reporting and logging

def test_python_errors_become_readable_tool_errors(allowed):
    tool = server.only_allowed_folder(server.h5ad_summary)
    with pytest.raises(ToolError, match="Column 'nope' not found"):
        tool("small.h5ad", column="nope")


def test_calls_are_logged_with_timing(allowed, caplog):
    tool = server.only_allowed_folder(server.vcf_stats)
    with caplog.at_level(logging.INFO, logger="biofile"):
        tool("test.vcf")
    messages = [r.getMessage() for r in caplog.records]
    assert any(m.startswith("call vcf_stats path=test.vcf") for m in messages)
    assert any(m.startswith("done vcf_stats in") for m in messages)


def test_refused_calls_are_logged_as_warnings(allowed, caplog):
    tool = server.only_allowed_folder(server.vcf_stats)
    with caplog.at_level(logging.INFO, logger="biofile"), pytest.raises(ToolError):
        tool("/etc/passwd")
    assert any(r.levelname == "WARNING" and "refused vcf_stats" in r.getMessage() for r in caplog.records)


def test_wrapper_keeps_name_and_docstring():
    tool = server.only_allowed_folder(server.bam_stats)
    assert tool.__name__ == "bam_stats"
    assert tool.__doc__ == server.bam_stats.__doc__


# End to end over stdio

def run_client(data_dir: Path, log_file: Path, calls: list[tuple[str, dict]]):
    params = StdioServerParameters(
        command=sys.executable,
        args=[str(PROJECT_DIR / "server.py")],
        cwd=str(PROJECT_DIR),
        env={**os.environ, "BIOFILE_ALLOWED_DIR": str(data_dir), "BIOFILE_LOG_FILE": str(log_file)},
    )

    async def main():
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                tools = (await session.list_tools()).tools
                results = [await session.call_tool(name, args) for name, args in calls]
                return tools, results

    return asyncio.run(main())


def test_mcp_client_end_to_end(data_dir, tmp_path):
    log_file = tmp_path / "biofile.log"
    tools, results = run_client(
        data_dir,
        log_file,
        [
            ("vcf_stats", {"path": "test.vcf"}),
            ("h5ad_summary", {"path": "small.h5ad", "column": "cell_type"}),
            ("vcf_stats", {"path": "/etc/passwd"}),
        ],
    )

    # Every tool is listed, read-only, and carries its docstring
    assert sorted(t.name for t in tools) == sorted(f.__name__ for f in server.TOOLS)
    assert len(tools) == 8
    assert all(t.annotations.read_only_hint for t in tools)
    assert all(t.description for t in tools)

    ok_vcf, ok_h5ad, refused = results
    assert not ok_vcf.is_error
    assert json.loads(ok_vcf.content[0].text)["ts_tv"] == 2.0
    assert json.loads(ok_h5ad.content[0].text)["cells_per_cell_type"]["B cell"] == 60

    assert refused.is_error
    assert "outside the allowed folder" in refused.content[0].text

    # If logs had gone to stdout the session above would have failed;
    # check they reached the log file instead
    log_text = log_file.read_text()
    assert "starting biofile server" in log_text
    assert "done vcf_stats" in log_text
    assert "refused vcf_stats" in log_text
