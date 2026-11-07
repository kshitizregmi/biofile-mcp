"""Take the web interface screenshots used in the format guides.

Opens each real data file in the running web interface, presses Summarize and
saves one picture per results tab into ``docs/screenshots/real/``.

Usage::

    streamlit run app.py --server.headless true --theme.base light
    python docs/make_screenshots.py           # in a second terminal

The example files must be in ``real_data/`` first; the download command is in
``docs/README.md``.

Needs Playwright with Chromium (``pip install playwright``, then
``playwright install chromium``). Pass the app address as the first argument
if it is not http://localhost:8501.
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

PROJECT_DIR = Path(__file__).resolve().parents[1]
DATA = PROJECT_DIR / "real_data"
OUT = PROJECT_DIR / "docs" / "screenshots" / "real"
URL = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8501"

# (file, region for BAM, column for h5ad, {tab name: picture name})
PAGES = [
    ("1_S103_L001_R1_001.fastq.gz", None, None,
     {"Quality": "fastq_quality", "Composition": "fastq_composition", "Read length": "fastq_read_length"}),
    ("test.paired_end.sorted.bam", "chr22:1500-5000", None,
     {"Overview": "bam_overview", "Coverage": "bam_coverage", "Mapping quality": "bam_mapq",
      "Insert size": "bam_insert_size", "Read length": "bam_read_length"}),
    ("NA12878_GIAB.chr22.vcf.gz", None, None,
     {"Overview": "vcf_overview", "Quality & depth": "vcf_quality_depth",
      "Indels": "vcf_indels", "Samples": "vcf_samples"}),
    ("pbmc3k.h5ad", None, "louvain",
     {"Groups": "h5ad_groups", "QC metrics": "h5ad_qc", "Embedding": "h5ad_embedding"}),
]


def wait_for_charts(page):
    """Wait until Streamlit has finished running and the charts are drawn."""
    page.wait_for_selector('[data-testid="stStatusWidget"]', state="detached", timeout=120_000)
    time.sleep(2)


OUT.mkdir(parents=True, exist_ok=True)

with sync_playwright() as p:
    browser = p.chromium.launch()
    # A tall window, so each tab fits on screen without scrolling
    page = browser.new_page(viewport={"width": 1400, "height": 3000}, device_scale_factor=1)

    for name, region, column, tabs in PAGES:
        page.goto(URL)
        page.get_by_text("Path on this computer").click()
        path_box = page.get_by_label("Full path")
        path_box.fill(str(DATA / name))
        path_box.press("Enter")
        wait_for_charts(page)

        if region:
            region_box = page.get_by_label("Region (optional)")
            region_box.fill(region)
            region_box.press("Enter")
            wait_for_charts(page)
        if column:
            page.get_by_label("Group cells by").click()
            page.get_by_role("option", name=column, exact=True).click()
            wait_for_charts(page)

        page.get_by_role("button", name="Summarize").click()
        wait_for_charts(page)

        for i, (tab, picture) in enumerate(tabs.items()):
            page.get_by_role("tab", name=tab, exact=True).click()
            wait_for_charts(page)
            panel = page.locator('[role="tabpanel"]:visible').bounding_box()
            # The first picture of each file also shows the file name, the
            # options and the headline numbers; the others show the tab only
            top = page.locator("h1").first.bounding_box()["y"] - 16 if i == 0 else panel["y"]
            clip = {"x": panel["x"] - 24, "y": top,
                    "width": panel["width"] + 48, "height": panel["y"] + panel["height"] + 16 - top}
            page.screenshot(path=str(OUT / f"{picture}.png"), clip=clip)
            print("saved", picture)

    browser.close()
