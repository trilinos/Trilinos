#!/bin/bash

exe() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed.";
    exit 1;
  fi
}

exe module load sierra-devel

SIERRA_CODE=$(git rev-parse --show-toplevel)
if [ $? -ne 0 ] ; then
  echo "Stopping due to error finding top level of sierra code repo.";
  exit 1;
fi

echo "SIERRA_CODE: ${SIERRA_CODE}"

exe cd $SIERRA_CODE
LOG_DIR=${SIERRA_CODE}/coupling_testing_logs
BIN_DIR_CURRENT=${SIERRA_CODE}/bin_current
BIN_DIR_PREVIOUS=${SIERRA_CODE}/bin_previous
BIN_DIR_TEMP=${SIERRA_CODE}/bin_temp
exe 'mkdir -p $LOG_DIR'
exe 'mkdir -p $BIN_DIR_CURRENT'
exe 'mkdir -p $BIN_DIR_PREVIOUS'
exe 'mkdir -p $BIN_DIR_TEMP'

branch_name="$(git symbolic-ref -q --short HEAD)" ||
branch_name="(unnamed branch)"
echo "$SIERRA_CODE is currently on branch $branch_name"

exe 'repo sync >& $LOG_DIR/repo_sync.log'
exe 'git rebase master'

exe 'assign -p stk_integration_tests*mock >& $LOG_DIR/assign_tests.log'
exe 'bake --bin-dir $BIN_DIR_CURRENT --targets-from-tests=assigned.tests >& $LOG_DIR/bake_current.log'

exe 'testrun --allow-multipliers time --bin-dir $BIN_DIR_CURRENT'

exe 'git checkout $(git rev-list -1 --before="1 day ago" master)'

exe_with_cleanup() {
  stdbuf -o0 -e0 echo "% $@" ;
  eval "$@" ;
  if [ $? -ne 0 ] ; then
    echo "'$@' failed, cleaning up.";
    exe 'rm -rf $BIN_DIR_PREVIOUS $BIN_DIR_CURRENT $BIN_DIR_TEMP'
    exe 'git checkout $branch_name'
    echo "Tests failed. Cleaned up bin directories and checked out branch $branch_name."
    exit 1;
  fi
}

exe_with_cleanup 'bake --bin-dir $BIN_DIR_PREVIOUS --targets-from-tests=assigned.tests >& $LOG_DIR/bake_previous.log'

exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_aria $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_salinas $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_fuego $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_sparc $BIN_DIR_TEMP'

echo 'Testing couplings with mock_aria and mock_salinas from today, and mock_fuego and mock_sparc from yesterday.'

exe_with_cleanup 'testrun --allow-multipliers time --bin-dir $BIN_DIR_TEMP --save-all-results'

exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_aria $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_salinas $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_fuego $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_sparc $BIN_DIR_TEMP'

echo 'Testing couplings with mock_fuego and mock_sparc from today, and mock_aria and mock_salinas from yesterday.'

exe_with_cleanup 'testrun --allow-multipliers time --bin-dir $BIN_DIR_TEMP --save-all-results'

exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_aria $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_PREVIOUS/mock_fuego $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_salinas $BIN_DIR_TEMP'
exe_with_cleanup '"cp" $BIN_DIR_CURRENT/mock_sparc $BIN_DIR_TEMP'

echo 'Testing couplings with mock_aria and mock_fuego from today, and mock_salinas and mock_sparc from yesterday.'

exe_with_cleanup 'testrun --allow-multipliers time --bin-dir $BIN_DIR_TEMP --save-all-results'

exe 'rm -rf $BIN_DIR_PREVIOUS $BIN_DIR_CURRENT $BIN_DIR_TEMP'
exe 'git checkout $branch_name'

echo "Tests passed. Cleaned up bin directories and checked out branch $branch_name."
