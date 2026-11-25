#!/usr/bin/env bash
# Move every submodule to the release tag that modules.versions pins, and
# stage the moves. Prints what moved; writes changed=0|1 to $GITHUB_OUTPUT
# when run in Actions.
set -euo pipefail

CHANGED=0
while IFS='=' read -r module version; do
  module=$(echo "$module" | xargs)
  version=$(echo "$version" | xargs)
  [ -z "$module" ] && continue
  [[ "$module" == \#* ]] && continue

  key=$(git config -f .gitmodules --get-regexp '^submodule\..*\.url$' \
    | awk -v m="$module" '$2 ~ ("/"m"(\\.git)?$") { sub(/\.url$/, ".path", $1); print $1 }' \
    | head -n1)
  if [ -z "$key" ]; then
    echo "WARN no submodule found for ${module}, skipping"
    continue
  fi
  path=$(git config -f .gitmodules --get "$key")
  tag="v${version#v}"

  git -C "$path" fetch --quiet --depth 1 origin "refs/tags/${tag}:refs/tags/${tag}"
  want=$(git -C "$path" rev-parse "refs/tags/${tag}^{commit}")
  if [ "$(git -C "$path" rev-parse HEAD)" = "$want" ]; then
    continue
  fi
  git -C "$path" checkout --quiet "$want"
  git add "$path"
  CHANGED=1
  echo "Moved ${path} to ${tag} (${want:0:8})"
done < modules.versions

echo "changed=${CHANGED}" >> "${GITHUB_OUTPUT:-/dev/null}"
[ "$CHANGED" = 1 ] || echo "All submodules already match modules.versions."
