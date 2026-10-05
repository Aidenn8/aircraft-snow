"""Assemble site/index.html from site/template.html, the replay and page data, and the scripts.

usage: build_site.py [--out site/index.html]
"""
import argparse
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "site" / "index.html"))
    ap.add_argument("--hero", default="")
    a = ap.parse_args()
    site = ROOT / "site"
    html = (site / "template.html").read_text()
    replay = json.loads((site / "denver_replay.json").read_text())
    html = html.replace("{{REPLAY_JSON}}", json.dumps(replay, separators=(",", ":")).replace("</", "<\\/"))
    html = html.replace("{{REPLAY_JS}}", (site / "replay.js").read_text())
    page = json.loads((site / "page_data.json").read_text())
    html = html.replace("{{PAGE_JSON}}", json.dumps(page, separators=(",", ":")).replace("</", "<\\/"))
    html = html.replace("{{CHARTS_JS}}", (site / "charts.js").read_text())
    html = html.replace("{{HERO_SUMMARY}}", a.hero)
    pathlib.Path(a.out).write_text(html)
    print(a.out, f"{len(html) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
