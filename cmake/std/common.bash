#!/usr/bin/env bash
#
# Common bash functions that we use in several scripts.
#

# Set some color codes if we're on a terminal that supports colors
if test -t 1; then
    # see if it supports colors...
    ncolors=$(tput colors)
    if test -n "$ncolors" && test $ncolors -ge 8; then
        export bold="$(tput bold)"
        export underline="$(tput smul)"
        export standout="$(tput smso)"
        export standout_end="$(tput rmso)"
        export normal="$(tput sgr0)"
        export dim="$(tput dim)"
        export black="$(tput setaf 0)"
        export red="$(tput setaf 1)"
        export green="$(tput setaf 2)"
        export yellow="$(tput setaf 3)"
        export blue="$(tput setaf 4)"
        export magenta="$(tput setaf 5)"
        export cyan="$(tput setaf 6)"
        export white="$(tput setaf 7)"
    fi
fi



# message_std
#
# param1: prefix
# param2: message text
function message_std
{
    local prefix=${1}
    local message=${2}
    echo -e "${prefix}${message}"
}



# print_centered_text
#
# Prints out a centered text string with endcaps
#
# param1: width
# param2: endcaps
# param3: text to print
function print_centered_text()
{
    local width=${1:?}
    local endcap=${2:?}
    local text=${3:?}
    local textsize=${#text}
    local capsize=${#endcap}
    local span1=$((($width + $textsize - $capsize * 2)/2))
    local span2=$(($width - $span1 - $capsize * 2))
    printf "%s%${span1}s%${span2}s%s\n" "${endcap}" "${text}" "" "${endcap}"
}



# print_banner()
#
# Prints out a banner block with date/time stamp.
#
# param1: banner text to print
function print_banner()
{
    local banner_text=${1:?}
    local textdate=$(date +"%Y-%m-%d %H:%M:%S")
    local width=60
    echo -e ""
    echo -e "+----------------------------------------------------------+"
    print_centered_text ${width} "|" "${banner_text}"
    print_centered_text ${width} "|" " "
    print_centered_text ${width} "|" "${textdate}"
    echo -e "+----------------------------------------------------------+"
}



# print_banner_2lines()
#
# Prints out a two line banner plus a date/time stamp.
# param1: banner text line 1
# param2: banner text line 2
function print_banner_2lines()
{
    local banner_text_line1=${1:?}
    local banner_text_line2=${2:?}
    local textdate=$(date +"%Y-%m-%d %H:%M:%S")
    local width=60
    echo -e ""
    echo -e "+----------------------------------------------------------+"
    print_centered_text ${width} "|" "${banner_text_line1}"
    print_centered_text ${width} "|" "${banner_text_line2}"
    print_centered_text ${width} "|" " "
    print_centered_text ${width} "|" "${textdate}"
    echo -e "+----------------------------------------------------------+"
}



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

