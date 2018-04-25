#!/bin/bash
user="GITHUB_USER_NAME"
fork="trilinos"
repo="Trilinos"
mainBranch="develop"

approval_string="Approved for push"
tokenfile=~/.githubOAuth/token

TMPFILE=/tmp/.ac$$

usage="Usage: isawesome.sh <PR>, where <PR> is a valid pull request number for $fork/$repo, "

# Ensure that we were only passed a single command line argument.
if [[ $# -ne 1 ]]; then
  echo $usage
  exit 1
fi

# Ensure that the argument was an integer.
re='^[0-9]+$'
if ! [[ $1 =~ $re ]]; then
  echo $usage
  echo "$1 must be an integer."
  exit 2
fi

# Get the details from that pull request.
pr=$1
token=$(cat $tokenfile)
h="'Authorization: token $token'"
curl -i -H $h https://api.github.com/repos/$fork/$repo/pulls/$pr >$TMPFILE 2> /dev/null

# Ensure the pull request number was valid.
if grep "message.*Not Found" $TMPFILE > /dev/null; then
  echo $usage
  echo "It looks like pull request $pr does not exist in $fork/$repo."
  exit 3
fi

# Get the SHA of the PR
SHA=`grep -A 3 '"head"' $TMPFILE | grep sha |cut -f2 -d: | sed 's/[", ]//g'`

# Make sure everyone knows this PR is awesome
CMD=$(echo curl -i -H \'Authorization: token $token\' -d \'{\"state\": \"success\" , \"description\": \"This code is awesome\", \"context\": \"awesome\"}\' https://api.github.com/repos/trilinos/Trilinos/statuses/$SHA)
eval $CMD > $TMPFILE 2> /dev/null

if grep state $TMPFILE > /dev/null; then
    echo "Auto-awesome succeeded"
else
    echo "Auto-awesome FAILED"; exit 3
fi

rm -f $TMPFILE