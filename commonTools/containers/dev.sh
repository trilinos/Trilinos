#!/bin/bash

# This script tries to install the Python dependencies before running dev-container
# This avoids having to provide instructions for installing dependencies.

SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
TRILINOS_SOURCE_DIR=${SCRIPT_DIR}/../..

if [[ ! -d ${SCRIPT_DIR}/.dev_venv || ! -f ${SCRIPT_DIR}/.dev_venv/succesful_install ]] ; then
    rm -rf ${SCRIPT_DIR}/.dev_venv
    echo "Setting up a virtual environment in ${SCRIPT_DIR}/.dev_venv"
    echo ""
    python3 -m venv ${SCRIPT_DIR}/.dev_venv
    source ${SCRIPT_DIR}/.dev_venv/bin/activate
    python3 -m pip install ${SCRIPT_DIR} && touch ${SCRIPT_DIR}/.dev_venv/succesful_install
else
    source ${SCRIPT_DIR}/.dev_venv/bin/activate
    python3 -m pip install ${SCRIPT_DIR} && touch ${SCRIPT_DIR}/.dev_venv/succesful_install
fi

if [[ -f ${SCRIPT_DIR}/.dev_venv/succesful_install ]] ; then
    # For example: ./dev.sh --container-engine podman
    dev-container "$@"

    deactivate
else
    deactivate

    echo ""
    echo "The installation of Python dependencies failed."
    echo "One possible reason is that the used Python is too old. This script needs Python 3.9 or newer. This system has"
    python3 --version

    echo ""
    echo "Another potential problem could be misconfigured certificates."
    echo "You could try to set PIP_CERT to the location of your certificate bundle"
    if [ -f /etc/ssl/certs/ca-bundle.crt ] ; then
        echo "Try"
        echo " export PIP_CERT=/etc/ssl/certs/ca-bundle.crt"
    fi
fi
