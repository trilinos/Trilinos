#!/bin/bash

# Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
# Solutions of Sandia, LLC (NTESS). Under the terms of Contract
# DE-NA0003525 with NTESS, the U.S. Government retains certain rights
# in this software.
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
#     * Neither the name of NTESS nor the names of its contributors
#       may be used to endorse or promote products derived from this
#       software without specific prior written permission.
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
STK_HEADER_FILE=${1:-stk_version.hpp}

# The current working directory might not be the source tree, so go there
# so that the git commands will work (if the source tree is a git repository.)
# Note that it seems a bad assumption that this script will also be in 
# the source tree, but I'm not going to change that now.
cd ${DIR}

REMOTE_REPO="git@cee-gitlab.sandia.gov:1540-compsim/sierra/base"

# Don't prepend DIR here, since we are now in that directory.
test -r $STK_HEADER_FILE &&
    CURRENT_VERSION=$(cat $STK_HEADER_FILE)

# If this is a git repo get its version. Otherwise we'll try other
# approaches to getting the version below, depending on whether this
# is an internal or external build.
NEW_VERSION=$(git describe --long --abbrev=8 --match=[0-9]*.[0-9]*.[0-9]* HEAD 2>/dev/null)

# If we do not have access to the REMOTE_REPO assume this is an external customer build.
if ! git ls-remote --exit-code --tags ${REMOTE_REPO} >/dev/null 2>&1
then
    # External user builds.
    # This could still be a clone of our repo (ie, Goodyear), in which case
    # NEW_VERSION will be set. If not try some fallbacks.
    test -z "${NEW_VERSION}" && NEW_VERSION=${CURRENT_VERSION##* }
    test -z "${NEW_VERSION}" && NEW_VERSION=unset
else
    # Internal user builds
    if [ -z "${NEW_VERSION}" ] ; then
        # Not in a local repo; try to use current header file.
        if [ -n "${CURRENT_VERSION:-}" ]
        then
            NEW_VERSION=${CURRENT_VERSION##* }
        else
            NEW_VERSION="no-version"
        fi
    else
        # In a local repo; check for changes.
        ROOT=$(git rev-parse --show-toplevel)
        test -z "$(git ls-files -m $ROOT)" ||
            NEW_VERSION="${NEW_VERSION}-modified"
    fi

    # Trim the 'v' from the beginning of the version (why is it there in the first place?)
    NEW_VERSION=$(expr "$NEW_VERSION" : v*'\(.*\)')
fi

# Don't touch the header if it already contains this version (to prevent unwanted recompilation.)
test "// $NEW_VERSION" = "$CURRENT_VERSION" ||
    echo "// $NEW_VERSION" > $STK_HEADER_FILE

echo -n "$NEW_VERSION"
