"""Assemble site/index.html from site/template.html, the replay and page data, and the scripts.

usage: build_site.py [--out site/index.html] [--data site/page_data.json]
                     [--gallery figures/census_gallery.png] [--hero TEXT]

The census gallery figure is embedded as a JPEG data URI (1600 px wide), so the page stays
a single file.
"""
import argparse
import base64
import io
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def gallery_uri(path: pathlib.Path, width: int = 1600) -> str:
    if not path.exists():
        print(f"warning: no gallery figure at {path}")
        return ""
    from PIL import Image
    im = Image.open(path).convert("RGB")
    if im.width > width:
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "JPEG", quality=84, optimize=True, progressive=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "site" / "index.html"))
    ap.add_argument("--hero", default="")
    ap.add_argument("--data", default=str(ROOT / "site" / "page_data.json"))
    ap.add_argument("--gallery", default=str(ROOT / "figures" / "census_gallery.png"))
    a = ap.parse_args()
    site = ROOT / "site"
    html = (site / "template.html").read_text()
    replay = json.loads((site / "denver_replay.json").read_text())
    html = html.replace("{{REPLAY_JSON}}", json.dumps(replay, separators=(",", ":")).replace("</", "<\\/"))
    html = html.replace("{{REPLAY_JS}}", (site / "replay.js").read_text())
    page = json.loads(pathlib.Path(a.data).read_text())
    html = html.replace("{{PAGE_JSON}}", json.dumps(page, separators=(",", ":")).replace("</", "<\\/"))
    html = html.replace("{{CHARTS_JS}}", (site / "charts.js").read_text())
    html = html.replace("{{HERO_SUMMARY}}", a.hero)
    html = html.replace("{{GALLERY_IMG}}", gallery_uri(pathlib.Path(a.gallery)))
    pathlib.Path(a.out).write_text(html)
    print(a.out, f"{len(html) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
