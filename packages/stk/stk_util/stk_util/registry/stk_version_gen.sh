#!/bin/bash
# Copyright (c) 2013, Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
# 
#     * Neither the name of Sandia Corporation nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
DIR=$(dirname $0)

# The current working directory might not be the source tree, so go there
# so that the git commands will work (if the source tree is a git repository.)
# Note that it seems a bad assumption that this script will also be in 
# the source tree, but I'm not going to change that now.
cd ${DIR}

# Don't prepend DIR here, since we are now in that directory.
STK_HEADER_FILE=stk_version.hpp
OVERRIDE_FILE="version"

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
    NEW_VERSION=$(git describe --abbrev=8 HEAD 2>/dev/null) &&
    case "$NEW_VERSION" in
        *$LF*) (exit 1) ;;
        [0-9]*)
            git update-index -q --refresh > /dev/null
            test -z "$(git diff-index --name-only HEAD --)" ||
            NEW_VERSION="${NEW_VERSION}-modified" ;;
        esac
then
    # PKN: leaving this here in case we want to muck
    # with the formatting later.
    #NEW_VERSION=$(echo "$NEW_VERSION" | sed -e 's/-/./g');
    NEW_VERSION=$NEW_VERSION

# If there's no git version, try the override file.
elif test -f "$OVERRIDE_FILE"
then
    NEW_VERSION=$(cat ${OVERRIDE_FILE}) ||
        { 
        echo >&2 "stk_version_gen.sh: Unable to read file $(pwd)/${OVERRIDE_FILE}" 
        exit 1
        }
    test "$NEW_VERSION" != 'ERROR' ||
        { 
        echo >&2 "stk_version_gen.sh: Invalid version 'ERROR' in $(pwd)/${OVERRIDE_FILE}. Deleting that file." 
        rm -f $(pwd)/${OVERRIDE_FILE}
        exit 2
        }
else
    # Give up (with an appropriate message.)
    type -p git >/dev/null && 
        echo >&2 "stk_version_gen.sh: No .git repository to determine git version and $(pwd)/${OVERRIDE_FILE} file not found." ||
        echo >&2 "stk_version_gen.sh: No git binary in PATH and $(pwd)/${OVERRIDE_FILE} file not found."
    exit 3
fi

# Trim the 'v' from the beginning of the version (why is it there in the first place?)
NEW_VERSION=$(expr "$NEW_VERSION" : v*'\(.*\)')

if test -r $STK_HEADER_FILE
then
    CURRENT_VERSION=$(cat $STK_HEADER_FILE)
else
    CURRENT_VERSION=unset
fi
# Don't touch the header if it already contains this version (to prevent unwanted recompilation.)
test "// $NEW_VERSION" = "$CURRENT_VERSION" ||
    echo "// $NEW_VERSION" > $STK_HEADER_FILE

# Echo out so build system can use as macro definition, but only if --record was not passed.
test -z "$RECORD" &&
    echo -n "$NEW_VERSION" ||
    echo -n "$NEW_VERSION" > $OVERRIDE_FILE
