#!/bin/sh
# Call this script from the 'bootstrap' script of each package
_COMMON_TOOL_BASE_DIR=$1;
cp $_COMMON_TOOL_BASE_DIR/buildTools/*.pl ./config/.
cp $_COMMON_TOOL_BASE_DIR/refactoring/string-replace.pl ./config/.
