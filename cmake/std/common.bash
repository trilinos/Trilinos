#!/usr/bin/env bash
#
# Common bash functions that we use in several scripts.
#


# Gets the current script name (full path + filename)
function get_scriptname() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    local script_file="${script_path}/${script_name:?}"
    echo "${script_file}"
}


# Gets the path to the current script (full path)
function get_scriptpath() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    echo "${script_path}"
}



# Get the md5sum of a filename.
# param1: filename
# returns: md5sum of the file.
function get_md5sum() {
    local filename=${1:?}
    local sig=$(md5sum ${filename:?} | cut -d' ' -f1)
    echo "${sig:?}"
}



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

