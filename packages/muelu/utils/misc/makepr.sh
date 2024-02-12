#!/bin/bash
fork="trilinos"
repo="Trilinos"
mainBranch="develop"

# Get the PR template from the environment, otherwise try to grab it from the current directory
if [[ $TRILINOS_PR_TEMPLATE ]]; then
    PR_TEMPLATE=$TRILINOS_PR_TEMPLATE
else
    T_PREFIX=`pwd | sed 's/\(.*\)Trilinos\(.*\)/\1/'`
    PR_TEMPLATE=${T_PREFIX}Trilinos/.github/PULL_REQUEST_TEMPLATE.md
fi


# Get the editor command from the environment (follows the svn variable default ordering)
if [[ $VISUAL ]]; then    
    EDITOR_CMD=$VISUAL
elif [[ $EDITOR ]]; then
    EDITOR_CMD=$EDITOR
fi

#TODO It looks like the github API allows one to get a single file.
#TODO Rather than use a local Trilinos source for the issue template,
#TODO one could also pull the template from github so that the script
#TODO always uses the most up-to-date version.

tokenfile=~/.githubOAuth/token

TMPFILE=$(mktemp /tmp/makepr.XXXXXX)

USAGE="Usage: `basename $0` [-hfrbleis] \"PR title\""
OPTDESCR="\n  -h     -- help\n  -f       -- fork [${fork}]\n  -r     -- repository [${repo}]\n  -b     -- branch [${mainBranch}]\n  -t
-- team [github package name for @mentions and labels, CASE-SENSITIVE]\n  -e
-- (r)eviewer [github handle]\n  -i
-- issue [generate a github issue with the following text]\n   -s
-- summary/description [first comment, ideally should reference github issue]\n   -F
-- force the PR generation even if the existing checkout isn't clean"
EXAMPLE_USAGE="Example: makepr.sh -t \"MueLu\" -t \"Xpetra\" -e \"jhux2\" -e \"csiefer2\" -s \"Fixes issue #666\" \"MueLu: implement nifty feature\""

LABELS="\"AT: AUTOMERGE\""
reviewers=""

# Parse command line options.
FORCE=0
while getopts hvFf:r:b:t:e:s:i: OPT; do
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
        F)
            FORCE=1
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
            unset EDITOR_CMD
            ;;
        i)
            ISSUE_TEXT=$OUTARG
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
    if [ $FORCE -eq 0 ]; then
        echo "This repo contains file diffs.  Cannot generate PR"
        exit -1
    else
        echo "This repo contains file diffs.  But you used -F.  Proceeding under the assumption of user competence."
    fi
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

# Get the oauth token
token=$(cat $tokenfile)
h="'Authorization: token $token'"

# Get title string
TITLE_STRING="$*"

# Generate an issue (and then immediately close it), if we want that
if [ -z ${ISSUE_TEXT+x} ]; then :
else
    ISSUE_BODY_TMPFILE=$(mktemp /tmp/issue_body.XXXXXX)
    # Open the issue
    echo "{\"title\": \"$TITLE_STRING\" , \"body\": \"$ISSUE_TEXT\"}"  > ${ISSUE_BODY_TMPFILE}

    CMD=$(echo curl -i -H $h -d @${ISSUE_BODY_TMPFILE} https://api.github.com/repos/$fork/$repo/issues)
    eval $CMD >$TMPFILE 2> $TMPFILE

    # Get the Issue number
    ISSUE_NUM=`grep number\": $TMPFILE | cut -f2 -d: | cut -f1 -d, | sed 's/ *//'`

    if grep 'number' $TMPFILE > /dev/null; then
        echo "Issue $ISSUE_NUM generated"
    else
        echo "Issue generation failed"
        echo "See $TMPFILE and ${ISSUE_BODY_TMPFILE}."
        exit 1
    fi

    # Close the issue
    echo "{\"state\": \"closed\"}"  > ${ISSUE_BODY_TMPFILE}
    CMD=$(echo curl -i -H $h -d @${ISSUE_BODY_TMPFILE} https://api.github.com/repos/$fork/$repo/issues/$ISSUE_NUM)
    eval $CMD >$TMPFILE 2> $TMPFILE

    if grep 'closed_by' $TMPFILE > /dev/null; then
        echo "Issue $ISSUE_NUM closed"
    else
        echo "Issue $ISSUE_NUM closing failed"
        echo "See $TMPFILE and ${ISSUE_BODY_TMPFILE}."
        exit 1
    fi

    PR_FIRST_COMMENT="$PR_FIRST_COMMENT\n See Issue #$ISSUE_NUM"

    rm -f $ISSUE_BODY_TMPFILE
fi

# Create the PR body from the Trilinos PR template.
# Insert the first comment from above.
# Insert carriage returns everywhere so that markdown renders it correctly.
# Remove the line with double quotes, as that screws up the JSON parsing.

# Generate a new pull request
PR_TEXT_TMPFILE=$(mktemp /tmp/pr_text.XXXXXX)
PR_BODY_TMPFILE=$(mktemp /tmp/pr_body.XXXXXX)

if [ -z ${EDITOR_CMD+x} ]; then
    PR_BODY=${PR_FIRST_COMMENT}

else 
    # TODO: Clean this up
    PR_BODY=`awk -v firstComment="${PR_FIRST_COMMENT}" -v teamMentions="${MENTIONS}" '/@trilinos/ {print teamMentions "\\\n"; next} /Please describe your changes in detail/ {print $0 "\\\n"; print firstComment "\\\n"; next} /^$/ {print; next} /PackageName:/ {print "the title with PackageName:.\\\n"; next} 1 {print $0 "\\\n"}' ${PR_TEMPLATE}`
    echo "$PR_BODY" > ${PR_TEXT_TMPFILE}
    $EDITOR_CMD $PR_TEXT_TMPFILE; 
    PR_BODY=`cat $PR_TEXT_TMPFILE`
fi

echo "{\"title\": \"$TITLE_STRING\" , \"head\": \"$REMOTE\" ,\"base\": \"$mainBranch\", \"body\": \"$PR_BODY\"}" > ${PR_BODY_TMPFILE}
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
rm -f $PR_TEXT_TMPFILE
