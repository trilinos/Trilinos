#!/bin/bash
#
# This script gets the merge commits from TriBITS on the mainline branch since
# the last snaphsot.

# Get the SHA1 of the commit from TriBITS for most recent snapshot
lastTribitsSnapshotSha1=$(git log --pretty='%s' -1 --grep="Automatic snapshot commit from tribits" -- cmake/tribits | awk '{print $NF}')
#echo "lastTribitsSnapshotSha1 = '${lastTribitsSnapshotSha1}'"

# Get the merge commit on the TriBITS mainline branch for ${lastTribitsSnapshotSha1}
tribitsMergeCommitForLastSnapshot=$(cd TriBITS ; git log --pretty="%H" --reverse master --not ${lastTribitsSnapshotSha1} | head -n 1)
#echo "tribitsMergeCommitForLastSnapshot = '${tribitsMergeCommitForLastSnapshot}'"
# NOTE: ${tribitsMergeCommitForLastSnapshot} was the merge commit that pulled
# the TriBITS topic branch containing ${lastTribitsSnapshotSha1} into the
# TriBITS 'master' branch.  Therefore, we only want to report merge commits
# (i.e. TriBITS PRs) that were merged after that.

# Get the merge commits since ${tribitsMergeCommitForLastSnapshot}
cd TriBITS 2>&1 >> /dev/null
git log --first-parent --pretty=format:'%Cgreen%h%Creset "%s"%nAuthor: %an <%ae>%nDate:   %ad (%cr)%n' master --not ${tribitsMergeCommitForLastSnapshot}
cd -  2>&1 >> /dev/null
