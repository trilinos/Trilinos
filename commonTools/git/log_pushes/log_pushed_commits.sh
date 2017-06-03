#!/bin/bash -e

#
# This script runs and records if new commits have been pushed to the Trilinos
# GitHub 'develop' branch since the last time it ran.  If it detects that
# there are new commits, then the top commit is logged by adding it to the top
# of the file ../TrilinosPushLog.txt.  This therefore records pushes to
# tracking branch (i.e. 'develop').
#
# This script gets run in a loop by loop_log_pushed_commits.sh
#

cd Trilinos/
git fetch
UPDATED_COMMITS=`git log --oneline HEAD..@{u}`
echo "Updated commits:"
echo "$UPDATED_COMMITS"
if [ "$UPDATED_COMMITS" != "" ] ; then
  echo "Has updated commits to log!"
  echo &> ../TrilinosLastPushEntry.txt
  date &>> ../TrilinosLastPushEntry.txt
  echo &>> ../TrilinosLastPushEntry.txt
  git log --pretty=fuller -1 @{u} &>> ../TrilinosLastPushEntry.txt
  echo &>> ../TrilinosLastPushEntry.txt
  echo "Commits pushed:" &>> ../TrilinosLastPushEntry.txt
  echo "$UPDATED_COMMITS" &>> ../TrilinosLastPushEntry.txt
  { cat ../TrilinosLastPushEntry.txt ; cat ../TrilinosPushLog.txt ; } > ../TrilinosPushLog.txt.new
  mv ../TrilinosPushLog.txt{.new,}
  git pull
fi
