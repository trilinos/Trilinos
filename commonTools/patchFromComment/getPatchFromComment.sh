#!/bin/bash

# This script tries to install the Python dependencies before running getPatchFromComment.py
# This avoids having to provide instructions for installing dependencies.

SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
TRILINOS_SOURCE_DIR=${SCRIPT_DIR}/../..

if [[ ! -d ${SCRIPT_DIR}/.getPatchFromComment_venv || ! -f ${SCRIPT_DIR}/.getPatchFromComment_venv/succesful_install ]] ; then
    rm -rf ${SCRIPT_DIR}/.getPatchFromComment_venv
    echo "Setting up a virtual environment in ${SCRIPT_DIR}/.getPatchFromComment_venv"
    echo ""
    python3 -m venv ${SCRIPT_DIR}/.getPatchFromComment_venv
    source ${SCRIPT_DIR}/.getPatchFromComment_venv/bin/activate
    python3 -m pip install -r ${SCRIPT_DIR}/requirements.txt && touch ${SCRIPT_DIR}/.getPatchFromComment_venv/succesful_install
else
    source ${SCRIPT_DIR}/.getPatchFromComment_venv/bin/activate
fi

if [[ -f ${SCRIPT_DIR}/.getPatchFromComment_venv/succesful_install ]] ; then
    python3 ${SCRIPT_DIR}/getPatchFromComment.py "$@"

    deactivate
else
    deactivate

    echo ""
    echo "The installation of Python dependencies failed."
    echo "One possible reason is that the used Python is too old. This script needs Python 3.9 or newer. This system has"
    python3 --version
fi
