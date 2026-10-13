"""Assemble docs/index.html (the project page, served by GitHub Pages) from site/template.html, the replay and page data, and the scripts.

usage: build_site.py [--out docs/index.html] [--data site/page_data.json]
                     [--panels figures/gallery_panels] [--hero TEXT]

The census gallery panels (fig_gallery.py --panels) are embedded as PNG data URIs with their
titles set in HTML, so the page stays a single file and the grid reflows on narrow screens.
"""
import argparse
import base64
import html as htmllib
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
HERO = ("Sixteen flights were matched to streaks along their wind-drifted paths; the same flights rotated to other "
        "bearings around the radar produce about four such matches. Run unattended at eight airports over three "
        "winters, the same test found such streaks on 58 days where 29 would be false alarms.")


def gallery_panels(folder: pathlib.Path) -> str:
    meta_path = folder / "panels.json"
    if not meta_path.exists():
        print(f"warning: no gallery panels in {folder}")
        return ""
    esc = htmllib.escape
    out = []
    for m in json.loads(meta_path.read_text()):
        uri = "data:image/png;base64," + base64.b64encode((folder / m["file"]).read_bytes()).decode()
        alt = (f"Radar reflectivity near {m['place']} on {m['time']}: a streak of snow along the wind-carried "
               f"path of {m['callsign']}")
        out.append(
            f'        <div class="panel"><img src="{uri}" width="660" height="660" loading="lazy" alt="{esc(alt)}">\n'
            f'          <p class="pt">{esc(m["place"])}<span>{esc(m["site"])}</span></p>\n'
            f'          <p class="pm">{esc(m["time"])}<br><span class="mono">{esc(m["callsign"])}</span> · '
            f'{esc(m["type"])}<br>{m["alt_ft"]:,.0f} ft · {m["temp_c"]}&nbsp;°C · {m["age_min"]} min after passing</p></div>')
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "index.html"))
    ap.add_argument("--hero", default=HERO)
    ap.add_argument("--data", default=str(ROOT / "site" / "page_data.json"))
    ap.add_argument("--panels", default=str(ROOT / "figures" / "gallery_panels"))
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
    html = html.replace("{{GALLERY_PANELS}}", gallery_panels(pathlib.Path(a.panels)))
    pathlib.Path(a.out).write_text(html)
    print(a.out, f"{len(html) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
