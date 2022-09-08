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
# Print a message with an optional prefix
#
# param1: prefix, use empty string for no prefix (can be empty)
# param2: message text (can be empty)
message_std() {
    local prefix=${1}
    local message=${2}
    echo -e "${prefix}${message}" 2>&1
}


# executable_exists()
# Determines if a file exists and is executable.
#
# param1: executable (with path if necessary)
function executable_exists()
{
    local cmd=${1:?}
    local output=1
    if [ ! command -v ${cmd:?} &> /dev/null ]; then
        output=0
    fi
    echo ${output:?}
}


# execute_command()
# Executes a command with some extra checking plus formatting and logging
#
# param1: command to execute
function execute_command()
{
    local command=${1:?}
    message_std "PRDriver> " "${magenta}${command:?}${normal}"

    local is_executable=$(executable_exists ${command:?})

    if [[ "${is_executable}" == "1" ]]; then
        local _start=`date +%s`
        eval ${command:?}
        local err=$?
        local _stop=`date +%s`
        local _runtime=$((_stop-_start))

        if [ $err -ne 0 ]; then
            message_std "PRDriver> " "${red}FAILED ${normal}(${_runtime} s)"
        else
            message_std "PRDriver> " "${green}OK ${normal}(${_runtime} s)"
        fi
    else
        message_std "PRDriver> " "${red}ERROR: command '${command:?}' is not executable"
        message_std "PRDriver> " "${red}FAILED${normal}"
    fi
    return $err
}


# execute_command_checked()
# Executes a command with extra checking with some formatting, etc.
# Will call exit() if the command fails.
#
# param1: command to execute
function execute_command_checked()
{
    local command=${1:?}
    message_std "PRDriver> " "${magenta}${command:?}${normal}"

    local is_executable=$(executable_exists ${command:?})

    if [[ "${is_executable}" == "1" ]]; then
        eval ${command:?}
        local err=$?
        if [ $err -ne 0 ]; then
            message_std "PRDriver> " "Command failed with status: ${err}"
            message_std "PRDriver> " ""
            exit $err
        else
            message_std "PRDriver> " "${green}OK${normal}"
        fi
    else
        message_std "PRDriver> " "${red}ERROR: command '${command:?}' is not executable"
        message_std "PRDriver> " "${red}FAILED${normal}"
        message_std "PRDriver> " ""
        exit 32
    fi
}



# print_centered_text()
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
    printf "%s%${span1}s%${span2}s%s\n" "${endcap}" "${text}" "" "${endcap}" 2>&1
}



# print_banner()
#dd
# Prints out a banner block with date/time stamp.
#
# param1: banner text to print
function print_banner()
{
    local banner_text=${1:?}
    local textdate=$(date +"%Y-%m-%d %H:%M:%S")
    local width=60
    echo -e "" 2>&1
    echo -e "+----------------------------------------------------------+" 2>&1
    print_centered_text ${width} "|" "${banner_text}"
    print_centered_text ${width} "|" " "
    print_centered_text ${width} "|" "${textdate}"
    echo -e "+----------------------------------------------------------+" 2>&1
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
    echo -e "" 2>&1
    echo -e "+----------------------------------------------------------+" 2>&1
    print_centered_text ${width} "|" "${banner_text_line1}"
    print_centered_text ${width} "|" "${banner_text_line2}"
    print_centered_text ${width} "|" " "
    print_centered_text ${width} "|" "${textdate}"
    echo -e "+----------------------------------------------------------+" 2>&1
}


# envvar_append_or_create()
# Append a value to an envvar if it exists or create the envvar if it does not.
#
#  $1 = envvar name
#  $2 = string to append
function envvar_append_or_create() {
    message_std "PRDriver> " "ENVVAR: Append '${2}' to '${1}'"
    # envvar $1 is not set
    if [[ ! -n "${!1+1}" ]]; then
        message_std "PRDriver> " "export ${1}=${2}"
        export ${1}="${2}"
    else
        message_std "PRDriver> " "export ${1}=${!1}:${2}"
        export ${1}="${!1}:${2}"
    fi
}


# envvar_prepend_or_create()
# Prepend a value to an envvar if it exists or create the envvar if it doesn't.
#
#  $1 = envvar name
#  $2 = string to prepend
function envvar_prepend_or_create() {
    message_std "PRDriver> " "ENVVAR: Prepend '${2}' to '${1}'"
    # envvar $1 is not set
    if [[ ! -n "${!1+1}" ]]; then
        message_std "PRDriver> " "export ${1}=${2}"
        export ${1}="${2}"
    else
        message_std "PRDriver> " "export ${1}=${2}:${!1}"
        export ${1}="${2}:${!1}"
    fi
}


# envvar_set_or_create()
# Set an envvar.
#
#  $1 = envvar name
#  $2 = string to prepend
function envvar_set_or_create() {
    message_std "PRDriver> " "ENVVAR: Set ${1} <- '${2}'"
    export ${1}="${2}"
}

# get_scriptname()
# Gets the current script name (full path + filename)
function get_scriptname() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    local script_file="${script_path}/${script_name:?}"
    echo "${script_file}"
}


# get_scriptpath()
# Gets the path to the current script (full path)
function get_scriptpath() {
    # Get the full path to the current script
    local script_name=`basename $0`
    local script_path=$(dirname $(readlink -f $0))
    echo "${script_path}"
}


# get_md5sum()
# Get the md5sum of a filename.
#
# param1: filename
#
# returns: md5sum of the file.
function get_md5sum() {
    local filename=${1:?}
    local sig=$(md5sum ${filename:?} | cut -d' ' -f1)
    echo "${sig:?}"
}


#
# Get pip
# - @param1 python_exe - the python executable to install PIP for
#function get_pip() {
#    local python_exe=${1:?}
#
#    echo -e "--- get_pip():"
#    echo -e "--- Python: ${python_exe:?}"
#
#    # fetch get-pip.py
#    echo -e "--- curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py"
#    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
#
#    get_pip_args=(
#        --user
#        --proxy="http://user:nopass@proxy.sandia.gov:80"
#        --no-setuptools
#        --no-wheel
#    )
#    echo -e ""
#    echo -e "--- ${python_exe:?} ./get-pip.py ${get_pip_args[@]}"
#    ${python_exe:?} ./get-pip.py ${get_pip_args[@]}
#}



# get_python_packages()
# Install required Python packages using pip.
# installs: configparser, mock, pytest, pytest-cov
#
# param1: pip_exe - the pip binary to use, i.e., pip3.
function get_python_packages() {
    local pip_exe=${1:?}

    message_std "PRDriver> " "pip executable: ${pip_exe:?}"

    pip_args=(
        configparser
        mock
        pytest
        pytest-cov
    )
    message_std "PRDriver> " ""
    # ${pip_exe:?} install --user ${pip_args[@]}
    ${pip_exe:?} install --user "${pip_args[@]}"
}


# executable_exists()
# - Determines if a file exists and is executable.
# param1: executable (with path if necessary)
function executable_exists()
{
    local cmd=${1:?}
    local output=1
    if [ ! command -v ${cmd:?} &> /dev/null ]; then
        output=0
    fi
    echo ${output:?}
}


# execute_command_checked
# - Attempt to run a command and call exit() if it fails.
# param1: command to execute
function execute_command_checked()
{
    local command=${1:?}
    echo -e "PRDriver> ${magenta}${command:?}${normal}"

    local is_executable=$(executable_exists ${command:?})

    if [[ "${is_executable}" == "1" ]]; then
        eval ${command:?}
        local err=$?
        if [ $err -ne 0 ]; then
            echo -e "PRDriver> Command failed with status: ${err}"
            echo -e " "
            exit $err
        else
            echo -e "PRDriver> ${green}OK${normal}"
        fi
    else
        echo -e "PRDriver> ${red}ERROR: command '${command:?}' is not executable"
        echo -e "PRDriver> ${red}FAILED${normal}"
        echo -e " "
        exit 32
    fi
}


# executable_exists()
# - Determines if a file exists and is executable.
# param1: executable (with path if necessary)
function executable_exists()
{
    local cmd=${1:?}
    local output=1
    if [ ! command -v ${cmd:?} &> /dev/null ]; then
        output=0
    fi
    echo ${output:?}
}


# execute_command_checked
# - Attempt to run a command and call exit() if it fails.
# param1: command to execute
function execute_command_checked()
{
    local command=${1:?}
    echo -e "PRDriver> ${magenta}${command:?}${normal}"

    local is_executable=$(executable_exists ${command:?})

    if [[ "${is_executable}" == "1" ]]; then
        eval ${command:?}
        local err=$?
        if [ $err -ne 0 ]; then
            echo -e "PRDriver> Command failed with status: ${err}"
            echo -e " "
            exit $err
        else
            echo -e "PRDriver> ${green}OK${normal}"
        fi
    else
        echo -e "PRDriver> ${red}ERROR: command '${command:?}' is not executable"
        echo -e "PRDriver> ${red}FAILED${normal}"
        echo -e " "
        exit 32
    fi
}

