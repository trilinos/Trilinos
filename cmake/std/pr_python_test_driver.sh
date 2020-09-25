#!/usr/bin/env bash

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


#
# Get pip
# - @param1 python_exe - the python executable to install PIP for
function get_pip() {
    local python_exe=${1:?}

    echo -e "--- Python: ${python_exe:?}"

    # fetch get-pip.py
    echo -e "--- curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py"
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py

    get_pip_args=(
        --user
        --proxy="http://wwwproxy.sandia.gov:80"
        --no-setuptools
        --no-wheel
    )
    echo -e ""
    echo -e "--- ${python_exe:?} ./get-pip.py ${get_pip_args[@]}"
    ${python_exe:?} ./get-pip.py ${get_pip_args[@]}
}


#
# Install Python pacakges using pip
#
# - @param1 pip_exe - the pip binary to use, i.e., pip3.
#
function get_python_packages() {
    local pip_exe=${1:?}

    echo -e "--- Pip   : ${pip_exe:?}"

    pip_args=(
        --use-feature=2020-resolver
        configparser
        mock
        pytest
        pytest-cov
    )
    echo -e "--- ${pip_exe:?} install --user ${pip_args[@]}"
    ${pip_exe:?} install --user ${pip_args[@]}
}



echo "================================================================================"
echo "="
echo "=                   P Y T H O N   T E S T   D R I V E R"
echo "="
echo "================================================================================"
echo ""
echo "----------------------------------------"
echo "-   P A C K A G E   S E T U P"
echo "----------------------------------------"
: "${1:?Param 1 ERROR: Python Executable not set or empty}"
: "${2:?Param 2 ERROR: Pip Executable not set or empty}"

python_exe=${1:?}
pip_exe=${2:?}

get_pip ${python_exe:?}
get_python_packages ${pip_exe:?}


echo -e ""
echo -e "----------------------------------------"
echo -e "-   E X E C U T E   T E S T S"
echo -e "----------------------------------------"
echo -e "--- ${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing"
${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing

err=$?
echo -e ""
echo -e "--- STATUS: ${err:?}"


# Some Cleanup
rm get-pip.py >& /dev/null
#rm -rf ~/.local &> /dev/null

# Exit with the status from the python testing step
exit ${err:?}

