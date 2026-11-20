# Settings and helpers shared by the krino Trilinos snapshot scripts (sourced, not run).
# See README_trilinos_snapshot.txt.

set -euo pipefail

TOOLS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# This Sierra checkout (the scripts live in $SIERRA/krino/tools) and, beside it, the Trilinos fork.
export SIERRA=${SIERRA:-$(cd "$TOOLS_DIR/../.." && pwd)}
export TRILINOS=${TRILINOS:-$(dirname "$SIERRA")/Trilinos}
# The Sierra commit to snapshot: merged master, unless overridden for a local check.
export SIERRA_REF=${SIERRA_REF:-origin/master}

export TRILINOS_UPSTREAM_URL=git@github.com:trilinos/Trilinos
export TRILINOS_BRANCH=develop
export SNAPSHOT_BRANCH=krino-snapshot

# NO_FETCH=1 skips fetching and NO_PUSH=1 skips pushing, for runs without GitHub or cee-gitlab access.
export NO_FETCH=${NO_FETCH:-0}
export NO_PUSH=${NO_PUSH:-0}
# ALLOW_OVERWRITE=1 lets a snapshot overwrite changes made to packages/krino in Trilinos since the last one.
export ALLOW_OVERWRITE=${ALLOW_OVERWRITE:-0}

run() {
  echo "% $*"
  "$@"
}

fail() {
  echo "Cannot proceed: $*" >&2
  exit 1
}

print_settings() {
  echo "SIERRA:     $SIERRA"
  echo "SIERRA_REF: $SIERRA_REF"
  echo "TRILINOS:   $TRILINOS"
}

# On machines with the Sierra modules, use their tools (e.g. a recent git); elsewhere use the system's.
load_sierra_devel_if_available() {
  if type module >/dev/null 2>&1 && module load sierra-devel >/dev/null 2>&1; then
    echo "Loaded module sierra-devel"
  fi
}

verify_clean_repo() {
  [[ -z $(git -C "$1" status --porcelain) ]] || fail "repository at $1 is not clean"
}

verify_upstream_remote() {
  local url
  url=$(git -C "$TRILINOS" remote get-url upstream 2>/dev/null) ||
    fail "$TRILINOS has no upstream remote; add it with: git -C $TRILINOS remote add upstream $TRILINOS_UPSTREAM_URL"
  [[ $url == "$TRILINOS_UPSTREAM_URL" ]] || fail "$TRILINOS remote upstream is $url, not $TRILINOS_UPSTREAM_URL"
}

# The GitHub owner of the fork, from the origin url (git@github.com:OWNER/Trilinos or https://github.com/OWNER/Trilinos).
fork_owner() {
  git -C "$TRILINOS" remote get-url origin | sed -E 's#^.*github\.com[:/]([^/]+)/.*$#\1#'
}
