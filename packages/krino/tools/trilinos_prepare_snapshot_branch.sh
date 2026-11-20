#!/usr/bin/env bash

# Starts a krino snapshot: makes the krino-snapshot branch of the Trilinos fork upstream develop.
# See README_trilinos_snapshot.txt.

source "$(dirname "${BASH_SOURCE[0]}")/trilinos_snapshot_common.sh"

print_settings
load_sierra_devel_if_available

cd "$TRILINOS"
verify_clean_repo "$TRILINOS"
verify_upstream_remote

[[ $NO_FETCH == 1 ]] || run git fetch --all --prune

# The previous snapshot was squash-merged, so its commits are not ancestors of develop: start the branch
# over from develop (not tracking it, so that a bare git push cannot go upstream) ...
run git checkout --no-track -B "$SNAPSHOT_BRANCH" "upstream/$TRILINOS_BRANCH"
git branch --unset-upstream 2>/dev/null || true # -B keeps tracking set up when the branch was made

# ... and replace the fork's branch, whether or not GitHub deleted it after the last merge, so that
# trilinos_create_snapshot_commit.sh can push normally.
if [[ $NO_PUSH == 1 ]]; then
  echo "NO_PUSH=1: not pushing $SNAPSHOT_BRANCH to origin"
else
  run git push --force-with-lease origin "$SNAPSHOT_BRANCH"
fi

echo "$SNAPSHOT_BRANCH is upstream/$TRILINOS_BRANCH ($(git log -1 --format='%h %ad' --date=short))."
echo "Next: $(dirname "${BASH_SOURCE[0]}")/trilinos_create_snapshot_commit.sh"
