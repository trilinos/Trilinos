#!/usr/bin/env bash
SCRIPTFILE=$(realpath $BASH_SOURCE)
SCRIPTPATH=$(dirname $SCRIPTFILE)
source ${SCRIPTPATH}/common.bash

#
# Prepare Environment
#

# The existence of PYTHONHOME can cause some odd problems
# if we're seeing errors that look like `No module named 'encodings'`
# then a usual fix is to try and unset PYTHONHOME.
# I think this is due to a conflict in modules loaded by SEMS and then
# our use of a specific Python that's different.
# - See: https://stackoverflow.com/q/38132755/2059999 for more information
unset PYTHONHOME


print_banner "P Y T H O N   T E S T   D R I V E R"
message_std "" ""
print_banner "P A C K A G E   S E T U P"
: "${1:?Param 1 ERROR: Python Executable not set or empty}"
: "${2:?Param 2 ERROR: Pip Executable not set or empty}"

python_exe=${1:?}
pip_exe=${2:?}

message_std "--- " "python: ${python_exe:?}"
message_std "--- " "pip   : ${pip_exe:?}"

# Remove get_pip usage -- we only really needed this for python2 since
# the SEMS python3 distros come with pip3 installed.
# get_pip ${python_exe:?}
get_python_packages "${python_exe} -m pip"

message_std "" ""
print_banner "E X E C U T E   T E S T S"
message_std "--- " "${python_exe:?} -m pytest --cov=. --cov-report term-missing"
${python_exe:?} -m pytest --cov=. --cov-report term-missing

err=$?
message_std "" ""
message_std "--- " "STATUS: ${err:?}"


# Some Cleanup
rm get-pip.py >& /dev/null
#rm -rf ~/.local &> /dev/null

# Exit with the status from the python testing step
exit ${err:?}
