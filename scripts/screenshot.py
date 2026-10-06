"""Screenshot a built page with headless Chrome, in light or dark mode, at a given width.

usage: screenshot.py PAGE.html OUT.png [--width 1280] [--height 1600] [--dark] [--frame]
The page is wrapped in the same minimal skeleton the artifact host adds. Headless Chrome
will not lay a window out narrower than about 500 px, so for phone widths use --frame: the
page is shown in an iframe of exactly --width inside a wider window.
"""
import argparse
import pathlib
import subprocess
import tempfile

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("page")
    ap.add_argument("out")
    ap.add_argument("--width", type=int, default=1280)
    ap.add_argument("--height", type=int, default=1600)
    ap.add_argument("--dark", action="store_true")
    ap.add_argument("--frame", action="store_true")
    a = ap.parse_args()
    body = pathlib.Path(a.page).read_text()
    theme = ' data-theme="dark"' if a.dark else ' data-theme="light"'
    html = (f'<!doctype html><html{theme}><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">'
            '<style>body{margin:0}img{max-width:100%}[hidden]{display:none!important}</style></head>'
            f'<body>{body}</body></html>')
    with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False) as f:
        f.write(html)
        tmp = f.name
    url, win = f"file://{tmp}", a.width
    if a.frame:
        outer = (f'<!doctype html><html><body style="margin:0;background:#888">'
                 f'<iframe src="file://{tmp}" style="border:0;display:block" width="{a.width}" '
                 f'height="{a.height}"></iframe></body></html>')
        with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False) as f:
            f.write(outer)
            url, win = f"file://{f.name}", max(a.width + 40, 640)
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--hide-scrollbars",
                    "--allow-file-access-from-files",
                    f"--window-size={win},{a.height}", "--virtual-time-budget=5000",
                    f"--screenshot={a.out}", url], check=True, capture_output=True)
    print(a.out)


if __name__ == "__main__":
    main()
