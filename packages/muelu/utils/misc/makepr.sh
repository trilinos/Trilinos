#!/bin/bash
fork="trilinos"
repo="Trilinos"
mainBranch="develop"

tokenfile=~/.githubOAuth/token

TMPFILE=/tmp/.ac$$

# Check for a message
if [ $# -eq 0 ]; then
    echo "A message must be included"
    exit -1
fi

# Make sure there are no diff'd files
git diff-index --quiet HEAD
if [ $? -ne 0 ]; then
    echo "This repo contains file diffs.  Cannot generate PR"
    exit -1
fi

# Make sure the local branch matches $mainBranch
CBRANCH=`git rev-parse --abbrev-ref HEAD`
if [ "$CBRANCH" != "$mainBranch" ]; then
    echo "The current branch is $CBRANCH not $mainBranch.  Cannot generate PR"
    exit -1
fi

# Get SHA1 of current HEAD (cut down to 7 characters)
SHA=`git rev-parse HEAD | cut -c1-7`
REMOTE=$USER-$SHA

# Push this branch to remote with a new name
git push origin $CBRANCH:$REMOTE

TITLE_STRING="Auto-PR for SHA $SHA"

# Generate a new pull request
MESSAGE="$*"
token=$(cat $tokenfile)
h="'Authorization: token $token'"
CMD=$(echo curl -i -H $h -d \'{\"title\": \"$TITLE_STRING\" , \"head\": \"$REMOTE\" ,\"base\": \"$mainBranch\", \"body\": \"$MESSAGE\"}\' https://api.github.com/repos/$fork/$repo/pulls)
eval $CMD >$TMPFILE 2> $TMPFILE

# Get the PR number
PRN=`grep number\": $TMPFILE | cut -f2 -d:`

if grep Created $TMPFILE > /dev/null; then
    echo "PR $PRN created successfully"
else
    echo "PR Generation failed"; 
    exit 1
fi


rm -f $TMPFILE