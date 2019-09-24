#!/bin/bash
fork="trilinos"
repo="Trilinos"
mainBranch="develop"
TRILINOS_SOURCE=/home/jhu/software/checkin/Trilinos

#TODO It looks like the github API allows one to get a single file.
#TODO Rather than use a local Trilinos source for the issue template,
#TODO one could also pull the template from github so that the script
#TODO always uses the most up-to-date version.

tokenfile=~/.githubOAuth/token

TMPFILE=$(mktemp /tmp/makepr.XXXXXX)

USAGE="Usage: `basename $0` [-hfrbles] \"PR title\""
OPTDESCR="\n  -h     -- help\n  -f       -- fork [${fork}]\n  -r     -- repository [${repo}]\n  -b     -- branch [${mainBranch}]\n  -t
-- team [github package name for @mentions and labels, CASE-SENSITIVE]\n  -e
-- (r)eviewer [github handle]\n  -s     -- summary/description [first comment, ideally should reference github issue]"
EXAMPLE_USAGE="Example: makepr.sh -t \"MueLu\" -t \"Xpetra\" -e \"jhux2\" -e \"csiefer2\" -s \"Fixes issue #666\" \"MueLu: implement nifty feature\""

LABELS="\"AT: AUTOMERGE\""
reviewers=""

# Parse command line options.
while getopts hvf:r:b:t:e:s: OPT; do
    case "$OPT" in
        h)
            echo -e $USAGE
            echo -e $OPTDESCR
            echo -e "\n$EXAMPLE_USAGE"
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
        t)
            teams=("${teams[@]}" $OPTARG)
            ;;
        e)
            if [ -z "$reviewers" ]; then
                reviewers="\"$OPTARG\""
            else
                reviewers="$reviewers,\"$OPTARG\""
            fi
            ;;
        s)
            PR_FIRST_COMMENT=$OPTARG
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

if [[ -z $PR_FIRST_COMMENT ]]; then
  PR_FIRST_COMMENT="Auto-PR for SHA $SHA"
fi

for tt in "${teams[@]}"; do
  lctt=`echo "$tt" | tr '[:upper:]' '[:lower:]'`
  MENTIONS="$MENTIONS @trilinos/$lctt"
  LABELS="$LABELS,\"pkg: $tt\""
done

# Create the PR body from the Trilinos PR template.
# Insert the first comment from above.
# Insert carriage returns everywhere so that markdown renders it correctly.
# Remove the line with double quotes, as that screws up the JSON parsing.
PR_BODY=`awk -v firstComment="${PR_FIRST_COMMENT}" -v teamMentions="${MENTIONS}" '/@trilinos/ {print teamMentions "\\\n"; next} /Please describe your changes in detail/ {print $0 "\\\n"; print firstComment "\\\n"; next} /^$/ {print; next} /PackageName:/ {print "the title with PackageName:.\\\n"; next} 1 {print $0 "\\\n"}' ${TRILINOS_SOURCE}/.github/PULL_REQUEST_TEMPLATE.md`

# Generate a new pull request
TITLE_STRING="$*"
PR_BODY_TMPFILE=$(mktemp /tmp/pr_body.XXXXXX)
echo "{\"title\": \"$TITLE_STRING\" , \"head\": \"$REMOTE\" ,\"base\": \"$mainBranch\", \"body\": \"$PR_BODY\"}" > ${PR_BODY_TMPFILE}
token=$(cat $tokenfile)
h="'Authorization: token $token'"

CMD=$(echo curl -i -H $h -d @${PR_BODY_TMPFILE} https://api.github.com/repos/$fork/$repo/pulls)

eval $CMD >$TMPFILE 2> $TMPFILE

# Get the PR number
PRN=`grep number\": $TMPFILE | cut -f2 -d: | cut -f1 -d, | sed 's/ *//'`

if grep Created $TMPFILE > /dev/null; then
    echo "PR $PRN created successfully"
else
    echo "PR Generation failed"; 
    echo "See $TMPFILE and ${PR_BODY_TMPFILE}."
    exit 1
fi

# Add labels
CMD=$(echo curl -i -H $h -d \'[$LABELS]\' https://api.github.com/repos/$fork/$repo/issues/$PRN/labels)
eval $CMD >$TMPFILE 2> $TMPFILE

if grep 'AT: AUTOMERGE' $TMPFILE > /dev/null; then
    echo "PR $PRN labeled as: $LABELS"
else
    echo "PR $PRN label failed: $LABELS"; 
    echo "See $TMPFILE and ${PR_BODY_TMPFILE}."
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
        echo "See $TMPFILE and ${PR_BODY_TMPFILE}."
        exit 1
    fi
fi


rm -f $TMPFILE
rm -f $PR_BODY_TMPFILE
