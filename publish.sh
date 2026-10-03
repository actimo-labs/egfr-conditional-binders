#!/usr/bin/env bash
# publish.sh — get the code on GitHub and the anchor map on the Hub.
# Run after epitope_prep.py has produced out/ from the real 6ARU.
set -euo pipefail

GH_USER="${GH_USER:?set GH_USER}"
HF_USER="${HF_USER:?set HF_USER}"
REPO="egfr-conditional-binders"
DATASET="egfr-domain3-anchor-map"

# ---------- sanity ----------
for f in out/residues.csv out/anchors.csv out/patches.json; do
  [[ -f "$f" ]] || { echo "missing $f — run epitope_prep.py first"; exit 1; }
done

python - <<'PY'
import json
d = json.load(open("out/patches.json"))
n = d["n_meeting_spec"]
print(f"spec-compliant patches: {n}")
if n == 0:
    print("  ! no patch meets the >=3 anchor spec. Publish anyway — a negative")
    print("    result on a real target is worth more than another pipeline README.")
else:
    p = next(p for p in d["patches"] if p["meets_spec_v2"])
    print("  best:", ", ".join(f"{a['resname']}{a['precursor_num']}({a['class']})"
                               for a in p["anchors"]),
          f"| span {p['max_span_A']} A")
PY

# ---------- code -> GitHub ----------
git init -q 2>/dev/null || true
cat > .gitignore <<'EOF'
out*/
*.pdb
*.pka
*.fasta
.venv/
__pycache__/
EOF
git add README.md epitope_prep.py .gitignore
git commit -qm "epitope prep and titratable-anchor mapper for EGFR domain III" || true
gh repo create "$GH_USER/$REPO" --public --source=. --push 2>/dev/null \
  || { git remote add origin "git@github.com:$GH_USER/$REPO.git" 2>/dev/null || true
       git push -u origin HEAD; }

# ---------- anchor map -> HF Hub ----------
python -m pip install -q huggingface_hub

# huggingface-cli was removed in huggingface_hub v1.0 (replaced by `hf`), so
# check auth through the library instead of guessing which CLI is installed.
python - <<'PY' || { echo "not logged in - run: hf auth login"; exit 1; }
import sys
from huggingface_hub import whoami
try:
    me = whoami()
except Exception as e:
    print(f"  auth check failed: {e}")
    sys.exit(1)
orgs = [o["name"] for o in me.get("orgs", [])]
print(f"  logged in as {me['name']}" + (f" | orgs: {', '.join(orgs)}" if orgs else ""))
PY

TMP=$(mktemp -d)
cp out/residues.csv out/anchors.csv out/patches.json "$TMP/"
cp DATASET_CARD.md "$TMP/README.md"
sed -i.bak "s|https://github.com/|https://github.com/$GH_USER/$REPO|" "$TMP/README.md"
rm -f "$TMP/README.md.bak"

# create_repo takes the full owner/name, so this works for a personal account
# and an organisation identically. huggingface-cli repo create does not.
HF_REPO="$HF_USER/$DATASET" SRC="$TMP" python - <<'PY'
import os
from huggingface_hub import create_repo, upload_folder
repo = os.environ["HF_REPO"]
create_repo(repo, repo_type="dataset", exist_ok=True)
upload_folder(
    repo_id=repo, repo_type="dataset", folder_path=os.environ["SRC"],
    commit_message="EGFR domain III titratable anchor map from 6ARU",
)
print("uploaded:", repo)
PY
rm -rf "$TMP"

echo
echo "code:    https://github.com/$GH_USER/$REPO"
echo "dataset: https://huggingface.co/datasets/$HF_USER/$DATASET"