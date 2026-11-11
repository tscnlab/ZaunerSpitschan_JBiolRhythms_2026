"""Check local HTML asset/link destinations without fetching external URLs."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit, unquote
import sys

root = Path(sys.argv[1] if len(sys.argv) > 1 else "report").resolve()
site_prefix = "/ZaunerSpitschan_JBiolRhythms_2026/"
class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links = []; self.ids = set()
    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a: self.ids.add(a["id"])
        for key in ("href", "src"):
            if key in a: self.links.append((tag, a[key]))

pages = {}
for p in root.rglob("*.html"):
    h = Links(); h.feed(p.read_text()); pages[p] = h
errors = []
for path, page in pages.items():
    for tag, link in page.links:
        u = urlsplit(link)
        if u.scheme or u.netloc or not link or link.startswith("javascript:"): continue
        target = unquote(u.path)
        if target.startswith(site_prefix): target = target[len(site_prefix):]; base = root
        elif target.startswith("/"): target = target.lstrip("/"); base = root
        else: base = path.parent
        dest = (base / target).resolve() if target else path
        if dest.is_dir(): dest /= "index.html"
        if not dest.exists(): errors.append(f"{path.relative_to(root)}: missing {link}")
        elif u.fragment and dest in pages and unquote(u.fragment) not in pages[dest].ids:
            errors.append(f"{path.relative_to(root)}: missing anchor {link}")
if errors:
    print("\n".join(sorted(set(errors)))); sys.exit(1)
print(f"Checked {len(pages)} HTML pages: local files and anchors resolve.")
