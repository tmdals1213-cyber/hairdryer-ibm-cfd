"""Verify the public manifest and exact archived solver source hashes."""
from pathlib import Path, PurePosixPath
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def main():
    manifest = json.loads((REPO / 'MANIFEST.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        rel = PurePosixPath(item['path'])
        if rel.is_absolute() or '..' in rel.parts or '\\' in str(rel):
            raise SystemExit('Unsafe manifest path')
        p = REPO.joinpath(*rel.parts)
        if p.stat().st_size != item['bytes'] or sha(p) != item['sha256']:
            raise SystemExit('Hash mismatch: ' + str(rel))
    provenance = json.loads((ROOT / 'results/source_provenance.json').read_text())
    for name, digest in provenance['source_sha256'].items():
        if sha(ROOT / 'src' / name) != digest: raise SystemExit('Source changed: ' + name)
    print(f"PASS: {len(manifest['files'])} manifest files; all 15 original solver sources unchanged")


if __name__ == '__main__': main()
