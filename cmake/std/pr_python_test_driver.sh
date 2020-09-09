#!/usr/bin/env bash

echo "================================================================================"
echo "="
echo "= P Y T H O N   P A C K A G E   S E T U P"
echo "="
echo "================================================================================"
python_exe=${1:?}
pip_exe=${2:?}

echo -e ">>> Python: ${python_exe:?}"
echo -e ">>> Pip   : ${pip_exe:?}"
echo -e ""

#python_exe="/projects/sierra/linux_rh7/install/Python/2.7.5/bin/python"
#pip_exe="${HOME:?}/.local/bin/pip2"
# pip_exe="${HOME:?}/pip"

# fetch get-pip.py
echo -e ">>> curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py"
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py

get_pip_args=(
    --user
    --proxy="http://wwwproxy.sandia.gov:80"
    --no-setuptools
    --no-wheel
)
echo -e ""
echo -e ">>> ${python_exe:?} ./get-pip.py ${get_pip_args[@]}"
${python_exe:?} ./get-pip.py ${get_pip_args[@]}


# prepend $HOME/.local to PATH
#export PATH="${pip_exe}:${PATH}"

pip_args=(
    --use-feature=2020-resolver
    configparser
    mock
    pytest
    pytest-cov
)
echo -e ">>> ${pip_exe:?} install --user ${pip_args[@]}"
${pip_exe:?} install --user ${pip_args[@]}


echo -e ">>> ${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing"
${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing

err=$?

# Some Cleanup
rm get-pip.py >& /dev/null
#rm -rf ~/.local &> /dev/null

# Exit with the status from the python testing step
exit ${err:?}

