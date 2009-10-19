#!/bin/bash

DIR=$(dirname $0)

# The current working directory might not be the source tree, so go there
# so that the git commands will work (if the source tree is a git repository.)
# Note that it seems a bad assumption that this script will also be in 
# the source tree, but I'm not going to change that now.
cd ${DIR}

# Don't prepend DIR here, since we are now in that directory.
SVF=sierra_version.hpp
OVERRIDE_FILE="version"

DEF_VER=vERROR

RECORD=
while test "$#" -ne 0
do
    case "$1" in
    --record)
        RECORD=1; shift ;;
        *)
        echo "Invalid argument: $1"
        exit 1
    esac
done

LF='
'

# First see if the source tree is a git project. If not then
# see if there is a version file (included in release tarballs).
# Finally, default.
if test $(git rev-parse --git-dir 2>/dev/null) &&
    VN=$(git describe --abbrev=8 HEAD 2>/dev/null) &&
    case "$VN" in
        *$LF*) (exit 1) ;;
        v[0-9]*)
            git update-index -q --refresh > /dev/null
            test -z "$(git diff-index --name-only HEAD --)" ||
            VN="$VN-modified" ;;
        esac
then
    # PKN: leaving this here in case we want to muck
    # with the formatting later.
    #VN=$(echo "$VN" | sed -e 's/-/./g');
    VN=$VN
elif test -f "$OVERRIDE_FILE"
then
    VN=$(cat ${OVERRIDE_FILE}) || VN="$DEF_VER"
else
    VN="$DEF_VER"
fi

VN=$(expr "$VN" : v*'\(.*\)')

if test -r $SVF
then
    VC=$(cat $SVF)
else
    VC=unset
fi
test "// $VN" = "$VC" ||
    echo "// $VN" > $SVF

# echo out so build system can use as macro definition
test -z "$RECORD" &&
    echo -n "$VN" ||
    echo -n "$VN" > $OVERRIDE_FILE
