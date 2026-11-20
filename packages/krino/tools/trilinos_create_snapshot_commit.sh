#!/usr/bin/env bash

# Commits Sierra krino (at SIERRA_REF, merged master by default) onto the krino-snapshot branch of the
# Trilinos fork and pushes it.  Rerun it to add commits to an open snapshot PR.
# See README_trilinos_snapshot.txt.

source "$(dirname "${BASH_SOURCE[0]}")/trilinos_snapshot_common.sh"

# Parts of Sierra krino that are not part of the Trilinos package: the Sierra framework bridge, the bjam
# and stand-alone (non-TriBITS) builds, and the editor and agent files.
SIERRA_ONLY_PATHS=(krino_sierra Jamfile SierraCmakeCode.cmake pykrino/SierraCmakeCode.cmake .clang-format AGENTS.md)

# The Sierra version of a commit, as stk_version_gen.sh gives it for a checkout of that commit.
sierra_version() {
  git -C "$SIERRA" describe --long --abbrev=8 --match='[0-9]*.[0-9]*.[0-9]*' "$1" | sed 's/^v*//'
}

# Reports how far Sierra stk is ahead of the stk that Trilinos has, since krino is built against it there.
report_stk_snapshot() {
  local subject stkCommit
  subject=$(git -C "$TRILINOS" log -1 --format=%s --grep='^STK: Snapshot' -- packages/stk)
  echo "Trilinos stk: $subject"
  stkCommit=$(sed -nE 's/.* from Sierra .*-g([0-9a-f]+).*/\1/p' <<<"$subject")
  if [[ -n $stkCommit ]] && git -C "$SIERRA" cat-file -e "$stkCommit^{commit}" 2>/dev/null; then
    echo "Sierra stk commits since then: $(git -C "$SIERRA" rev-list --count "$stkCommit..$sierraCommit" -- stk)" \
      "(if krino uses any of them, the Trilinos build fails until the next stk snapshot)"
  fi
}

# Puts the files git tracks in Sierra krino at the commit (not ignored build or editor files, e.g. .cache,
# compile_commands.json, __pycache__) in directory/krino.
extract_tracked_sierra_krino() {
  git -C "$SIERRA" archive --format=tar "$1" krino | tar -x -C "$2"
}

# Puts the Trilinos package of Sierra krino at the commit in directory/krino.
extract_sierra_krino() {
  extract_tracked_sierra_krino "$1" "$2"
  local path
  for path in "${SIERRA_ONLY_PATHS[@]}"; do
    rm -rf "$2/krino/$path"
  done
}

update_krino_package() {
  echo "% replace packages/krino with Sierra krino at $1"
  rm -rf packages/krino
  extract_sierra_krino "$1" packages
  run git add --all packages/krino
}

# The newest commit of packages/krino that is a snapshot, and the Sierra version of its last iteration (a
# squash-merged snapshot PR lists every iteration's subject in its message).
last_snapshot_commit() {
  git log -1 --format=%H -E --grep='Krino: Snapshot .* from Sierra ' -- packages/krino
}

last_snapshot_sierra_version() {
  git log -1 --format=%B "$1" | sed -nE 's/.*Krino: Snapshot .* from Sierra ([^ ]+).*/\1/p' | tail -1
}

# The Sierra commit of a version, if it is unmodified and in this Sierra clone.
sierra_commit_of_version() {
  local hash
  hash=$(sed -nE 's/^.*-g([0-9a-f]+)$/\1/p' <<<"$1")
  [[ -n $hash ]] && git -C "$SIERRA" rev-parse --quiet --verify "$hash^{commit}"
}

# Whether two files have the same contents, or are both missing.
is_same_file() {
  [[ ! -e $1 && ! -e $2 ]] || cmp -s "$1" "$2"
}

# The files (relative to krino) that differ between, or are in only one of, two krino directories.
differing_files() {
  local file
  { (cd "$1" && find . -type f -o -type l); (cd "$2" && find . -type f -o -type l); } | sed 's#^\./##' | sort -u |
    while read -r file; do
      is_same_file "$1/$file" "$2/$file" || echo "$file"
    done
}

# Stops (unless ALLOW_OVERWRITE=1) if the snapshot would overwrite changes made to packages/krino in Trilinos
# since the last snapshot: files that differ from that snapshot's Sierra krino, and from the new one.  The last
# snapshot is compared with all its tracked files, whatever was left out of it, so that a change to
# SIERRA_ONLY_PATHS is not taken for a change made in Trilinos: a file left out then is missing from Trilinos
# and from the new snapshot, and one left out now is the same in Trilinos as in the last snapshot.
check_for_overwritten_trilinos_changes() {
  local lastSnapshot lastVersion lastSierraCommit
  lastSnapshot=$(last_snapshot_commit)
  if [[ -z $lastSnapshot ]]; then
    report_possible_overwrites "no earlier snapshot of packages/krino was found" ""
    return
  fi
  lastVersion=$(last_snapshot_sierra_version "$lastSnapshot")
  echo "Last snapshot: $(git log -1 --format='%h %ad' --date=short "$lastSnapshot"), from Sierra $lastVersion"
  if ! lastSierraCommit=$(sierra_commit_of_version "$lastVersion"); then
    local commitsSince
    commitsSince=$(commits_of_krino_since "$lastSnapshot")
    [[ -z $commitsSince ]] ||
      report_possible_overwrites "Sierra $lastVersion cannot be rebuilt to compare, and packages/krino has commits since" \
        "$commitsSince"
    return
  fi

  local workDir changedInTrilinos overwritten file
  workDir=$(mktemp -d)
  mkdir "$workDir/last" "$workDir/new"
  extract_tracked_sierra_krino "$lastSierraCommit" "$workDir/last"
  extract_sierra_krino "$sierraCommit" "$workDir/new"
  changedInTrilinos=$(differing_files "$workDir/last/krino" packages/krino)
  overwritten=$(while read -r file; do
      [[ -z $file ]] || is_same_file "$workDir/new/krino/$file" "packages/krino/$file" || echo "  $file"
    done <<<"$changedInTrilinos")
  rm -rf "$workDir"

  [[ -z $overwritten ]] ||
    report_possible_overwrites "these files were changed in Trilinos since the last snapshot, and differ in Sierra" \
      "$overwritten"$'\n'"Commits of packages/krino since the last snapshot:"$'\n'"$(commits_of_krino_since "$lastSnapshot")"
}

commits_of_krino_since() {
  git log --format='  %h %ad %an: %s' --date=short "$1..HEAD" -- packages/krino
}

# Stops unless ALLOW_OVERWRITE=1.
report_possible_overwrites() {
  echo "WARNING: the snapshot may overwrite changes made in Trilinos: $1"
  [[ -z $2 ]] || echo "$2"
  [[ $ALLOW_OVERWRITE == 1 ]] ||
    fail "port the changes to Sierra (or decide to drop them), then rerun with ALLOW_OVERWRITE=1 to overwrite them"
  echo "ALLOW_OVERWRITE=1: overwriting them"
}

print_settings
load_sierra_devel_if_available

[[ $NO_FETCH == 1 ]] || run git -C "$SIERRA" fetch origin
sierraCommit=$(git -C "$SIERRA" rev-parse --verify "$SIERRA_REF^{commit}") || fail "no Sierra commit $SIERRA_REF"
if ! git -C "$SIERRA" merge-base --is-ancestor "$sierraCommit" origin/master; then
  echo "WARNING: $SIERRA_REF is not merged into Sierra origin/master; do not open a PR from this snapshot."
fi
sierraVersion=$(sierra_version "$sierraCommit")
echo "Sierra krino: $sierraVersion"

cd "$TRILINOS"
verify_clean_repo "$TRILINOS"
[[ $(git branch --show-current) == "$SNAPSHOT_BRANCH" ]] ||
  fail "$TRILINOS is not on $SNAPSHOT_BRANCH; run trilinos_prepare_snapshot_branch.sh first"
report_stk_snapshot
check_for_overwritten_trilinos_changes

update_krino_package "$sierraCommit"
if git diff --cached --quiet; then
  echo "Trilinos krino already matches Sierra $sierraVersion: nothing to commit."
  exit 0
fi
git diff --cached --shortstat

run git commit --quiet -s -m "Krino: Snapshot $(date +'%m-%d-%y %H:%M') from Sierra $sierraVersion"
git log -1 --format='Committed %h %s'

if [[ $NO_PUSH == 1 ]]; then
  echo "NO_PUSH=1: not pushing $SNAPSHOT_BRANCH to origin"
  exit 0
fi
run git push origin "$SNAPSHOT_BRANCH"
echo "If there is no PR yet, open one against trilinos/Trilinos $TRILINOS_BRANCH, titled with the commit subject:"
echo "  https://github.com/trilinos/Trilinos/compare/$TRILINOS_BRANCH...$(fork_owner):Trilinos:$SNAPSHOT_BRANCH"
