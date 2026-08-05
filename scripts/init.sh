#!/usr/bin/env bash
# Scaffold this template into a concrete nf-<name> pipeline.
# Usage: ./scripts/init.sh <name>

set -e

NAME="$1"
if [ -z "$NAME" ]; then
  echo "usage: $0 <name>"
  exit 1
fi

LOWER=$(echo "$NAME" | tr '[:upper:]' '[:lower:]')
UPPER=$(echo "$NAME" | tr '[:lower:]' '[:upper:]')

# Files containing the template tokens.
FILES="main.nf README.md nextflow.config conf/base.config conf/test.config"

# Replace __NAME__ and __NAME_UPPER__ in each file.
for f in $FILES; do
  sed -e "s/__NAME_UPPER__/$UPPER/g" -e "s/__NAME__/$LOWER/g" "$f" > "$f.tmp"
  mv "$f.tmp" "$f"
done

# Strip the template-usage block from README.md.
sed '/<!-- TEMPLATE:START -->/,/<!-- TEMPLATE:END -->/d' README.md > README.md.tmp
mv README.md.tmp README.md

# Remove template-only scaffolding from the new repo.
git rm -f --ignore-unmatch scripts/init.sh .github/workflows/init.yml
rmdir scripts 2>/dev/null || true

echo "Scaffolded nf-$LOWER."
