#!/bin/bash
user="GITHUB_USER_NAME"
fork="trilinos"
repo="Trilinos"
mainBranch="develop"

approval_string="Approved for push"
tokenfile=~/.githubOAuth/token

TMPFILE=/tmp/.ac$$

usage="Usage: autoretester.sh <PR>, where <PR> is a valid pull request number for $fork/$repo, "



function msg() {
    D=`date "+%D %H:%M:%S"`
    echo "[$D] $*"
}


function error_check() {
    if grep 'API rate limit exceeded' $TMPFILE > /dev/null; then
        msg "WARNING: API rate limit exceeded"
    fi
}



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


# Has it been merged yet?
# NOTE: This check is somewhat different from what the API suggests
curl -i -H "Accept: application/vnd.github.v3+json"  https://api.github.com/repos/$fork/$repo/pulls/$pr/merge >$TMPFILE 2> /dev/null
if grep 'HTTP/1.1 204' $TMPFILE > /dev/null; then 
    merged=1
else
    merged=0
fi


while [ $merged -eq 0 ]; do    
    msg "PR $pr has NOT been merged"

    # Check for AT: RETEST
    curl -i -H $h https://api.github.com/repos/$fork/$repo/issues/$pr/labels >$TMPFILE 2> /dev/null
    error_check


    if grep '"name"' $TMPFILE| grep 'AT: RETEST' > /dev/null; then
        msg "PR $pr does have AT: RETEST applied"

    else
        msg "PR $pr does NOT have AT: RETEST applied... adding"

        # Add AT: RETEST
        CMD=$(echo curl -i -X POST -H \'Authorization: token $token\' -d \'{\"labels\": [\"AT: RETEST\"]}\' https://api.github.com/repos/$fork/$repo/issues/$pr/labels)
        eval $CMD > $TMPFILE 2> /dev/null
    fi

    # Sleep for an hour
    sleep 1800
    
    # Check to see if we're merged
    # NOTE: This check is somewhat different from what the API suggests
    curl -i -H "Accept: application/vnd.github.v3+json" https://api.github.com/repos/$fork/$repo/pulls/$pr/merge >$TMPFILE 2> /dev/null
    if grep 'HTTP/1.1 204' $TMPFILE > /dev/null; then 
        merged=1
    fi
    
done


# Final cleanup
msg "PR $pr has been MERGED!"        


# Make sure everyone knows this PR is awesome
#CMD=$(echo curl -i -H \'Authorization: token $token\' -d \'{\"state\": \"success\" , \"description\": \"This code is awesome\", \"context\": \"awesome\"}\' https://api.github.com/repos/trilinos/Trilinos/statuses/$SHA)

#rm -f $TMPFILE
