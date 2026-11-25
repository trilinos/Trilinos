Syncing Sierra krino to Trilinos (a "krino snapshot")

The snapshot is a PR from the krino-snapshot branch of a GitHub fork of Trilinos
(drnobleabq/Trilinos) to the develop branch of trilinos/Trilinos, replacing packages/krino with
Sierra krino as merged to Sierra master.

ONE-TIME SETUP

Clone the fork beside the Sierra checkout, so that the scripts find it, and add trilinos/Trilinos as
the upstream remote:

  cd <directory containing the Sierra checkout>
  git clone git@github.com:drnobleabq/Trilinos
  git -C Trilinos remote add upstream git@github.com:trilinos/Trilinos

The fork's own master and develop are not used; there is no need to pull or update them.

EACH SNAPSHOT

These need GitHub (and cee-gitlab, to fetch Sierra), so run them where both are reachable.

1. Merge the krino changes to Sierra master.  The snapshot is of Sierra origin/master, fetched by
   the script; the Sierra working tree is not checked out or changed, so it can be on any branch,
   with local changes.

2. Start the snapshot branch from upstream develop:

     <Sierra>/krino/tools/trilinos_prepare_snapshot_branch.sh

   This fetches all remotes, makes krino-snapshot upstream/develop (discarding the previous
   snapshot's branch, since snapshot PRs are squash-merged) and force-pushes it to the fork.

3. Commit and push the snapshot:

     <Sierra>/krino/tools/trilinos_create_snapshot_commit.sh

   This replaces packages/krino with the files git tracks in Sierra krino at origin/master (not
   ignored build or editor files), except those that are not part of the Trilinos package:
   krino_sierra, Jamfile, SierraCmakeCode.cmake and pykrino/SierraCmakeCode.cmake (so the
   stand-alone, non-TriBITS branch of krino's CMakeLists.txt files cannot be used from Trilinos),
   .clang-format and AGENTS.md.  It commits, signed off, as
     Krino: Snapshot MM-DD-YY HH:MM from Sierra <version>
   and pushes.  It also reports how many Sierra stk commits Trilinos's stk lacks: krino is built
   against Trilinos's stk, so a snapshot that uses newer stk must wait for the next stk snapshot.

   Before replacing packages/krino, it checks that the snapshot does not overwrite changes made to
   it in Trilinos since the last snapshot.  It rebuilds the last snapshot from the Sierra commit
   named in its message, and lists each file that Trilinos has changed, added or deleted since,
   and that differs from Sierra now (changes Sierra already has are not listed), with the
   Trilinos commits of packages/krino since.  If the last snapshot cannot be rebuilt (a
   "-modified" version, or a commit not in this clone), it lists the commits instead.  It then
   stops, leaving Trilinos unchanged: port the changes to Sierra, or decide to drop them, and
   rerun with ALLOW_OVERWRITE=1.

4. Open the PR with the link the script prints: base trilinos/Trilinos develop (not master),
   titled with the commit subject.  If more changes are needed, merge them to Sierra master and
   rerun step 3; each rerun adds a commit to the open PR.

OPTIONS (environment variables)

  SIERRA      Sierra checkout (default: the one containing these scripts)
  TRILINOS    Trilinos fork checkout (default: Trilinos beside $SIERRA)
  SIERRA_REF  Sierra commit to snapshot (default origin/master).  Another ref is for local checks
              only: the script warns that it is not merged, and its PR should not be opened.
  NO_FETCH=1  do not fetch (Sierra in create, all remotes in prepare)
  NO_PUSH=1   do not push
  ALLOW_OVERWRITE=1  overwrite changes made to packages/krino in Trilinos since the last snapshot

For example, to make a local snapshot commit of an unmerged branch for a test build of Trilinos on a
machine without GitHub or cee-gitlab access (afterwards rerun step 2 before a real snapshot):

  NO_FETCH=1 NO_PUSH=1 <Sierra>/krino/tools/trilinos_prepare_snapshot_branch.sh
  NO_FETCH=1 NO_PUSH=1 SIERRA_REF=HEAD <Sierra>/krino/tools/trilinos_create_snapshot_commit.sh
