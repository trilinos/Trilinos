#!/bin/bash -e

#
# This script runs and records if new commits have been pushed to the Trilinos
# GitHub 'develop' branch since the last time it ran.  If it detects that
# there are new commits, then the top commit is logged and appended to the
# file ../TrilinosPushLog.txt.  This therefore records pushes to tracking
# branch (i.e. 'develop').
#
# This script gets run in a loop by loop_log_pushed_commits.sh
#

cd Trilinos/
git fetch
UPDATED_COMMITS=`git log --oneline HEAD..@{u}`
echo "Updated commits:"
echo $UPDATED_COMMITS
if [ "$UPDATED_COMMITS" != "" ] ; then
  echo "Has updated commits to log!"
  echo &>> ../TrilinosPushLog.txt
  date &>> ../TrilinosPushLog.txt
  echo &>> ../TrilinosPushLog.txt
  git log --pretty=fuller -1 @{u} &>> ../TrilinosPushLog.txt
  echo &>> ../TrilinosPushLog.txt
  echo "Commits pushed:" &>> ../TrilinosPushLog.txt
  echo "$UPDATED_COMMITS" &>> ../TrilinosPushLog.txt
  git pull
fi
