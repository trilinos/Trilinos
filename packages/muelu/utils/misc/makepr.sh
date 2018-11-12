#!/bin/bash
fork="trilinos"
repo="Trilinos"
mainBranch="develop"

tokenfile=~/.githubOAuth/token

TMPFILE=/tmp/.ac$$

USAGE="Usage: `basename $0` [-hfrbl] \"PR title\""
OPTDESCR="\n  -h     -- help\n  -f       -- fork [${fork}]\n  -r     -- repository [${repo}]\n  -b     -- branch [${mainBranch}\n  -l     -- label [label_name]\n  -e     -- (r)eviewer [reviewer_name]"

labels="\"AT: AUTOMERGE\""
reviewers=""

# Parse command line options.
while getopts hvf:r:b:l:e: OPT; do
    case "$OPT" in
        h)
            echo -e $USAGE
            echo -e $OPTDESCR
            exit 0
            ;;
        v)
            echo "`basename $0` version 0.1"
            exit 0
            ;;
        f)
            fork=$OPTARG
            ;;
        r)
            repo=$OPTARG
            ;;
        b)
            mainBranch=$OPTARG
            ;;
        l)
            labels="$labels,\"$OPTARG\""
            ;;
        e)
            if [ -z "$reviewers" ]; then
                reviewers="\"$OPTARG\""
            else
                reviewers="$reviewers,\"$OPTARG\""
            fi
            ;;
        \?)
            # getopts issues an error message
            echo $USAGE >&2
            exit 1
            ;;
    esac
done

# Remove the options we parsed above.
shift `expr $OPTIND - 1`

# We want at least one non-option argument.
# Remove this block if you don't need it.
if [ $# -eq 0 ]; then
    echo -e $USAGE >&2
    echo -e $OPTDESCR
    exit 1
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

MESSAGE_STRING="Auto-PR for SHA $SHA"

# Generate a new pull request
TITLE_STRING="$*"
token=$(cat $tokenfile)
h="'Authorization: token $token'"
CMD=$(echo curl -i -H $h -d \'{\"title\": \"$TITLE_STRING\" , \"head\": \"$REMOTE\" ,\"base\": \"$mainBranch\", \"body\": \"$MESSAGE_STRING\"}\' https://api.github.com/repos/$fork/$repo/pulls)
eval $CMD >$TMPFILE 2> $TMPFILE

# Get the PR number
PRN=`grep number\": $TMPFILE | cut -f2 -d: | cut -f1 -d, | sed 's/ *//'`

if grep Created $TMPFILE > /dev/null; then
    echo "PR $PRN created successfully"
else
    echo "PR Generation failed"; 
    exit 1
fi

# Add labels
CMD=$(echo curl -i -H $h -d \'[$labels]\' https://api.github.com/repos/$fork/$repo/issues/$PRN/labels)
eval $CMD >$TMPFILE 2> $TMPFILE

if grep 'AT: AUTOMERGE' $TMPFILE > /dev/null; then
    echo "PR $PRN labeled as: $labels"
else
    echo "PR $PRN label failed: $labels"; 
    cat $TMPFILE
    exit 1
fi

# Add reviewers
if [ -z "$reviewers" ]; then :
else
    CMD=$(echo curl -i -H $h -d \'{\"reviewers\": [$reviewers]}\' https://api.github.com/repos/$fork/$repo/pulls/$PRN/requested_reviewers)
    eval $CMD >$TMPFILE 2> $TMPFILE
    
    if grep 'review_comments_url' $TMPFILE > /dev/null; then
        echo "PR $PRN adding reviewers : $reviewers"
    else
        echo "PR $PRN adding reviewers failed: $reviewers"; 
        cat $TMPFILE
        exit 1
    fi
fi


rm -f $TMPFILE