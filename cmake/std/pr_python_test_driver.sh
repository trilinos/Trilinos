#!/usr/bin/env bash

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

echo -e ">>> Python: ${python_exe:?}"
echo -e ">>> Pip   : ${pip_exe:?}"
echo -e ""

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


pip_args=(
    --use-feature=2020-resolver
    configparser
    mock
    pytest
    pytest-cov
)
echo -e ">>> ${pip_exe:?} install --user ${pip_args[@]}"
${pip_exe:?} install --user ${pip_args[@]}

echo -e ""
echo -e "----------------------------------------"
echo -e "-   E X E C U T E   T E S T S"
echo -e "----------------------------------------"
echo -e ">>> ${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing"
${python_exe:?} -m pytest --cov=trilinosprhelpers --cov-report term-missing

err=$?
echo -e ""
echo -e ">>> STATUS: ${err:?}"


# Some Cleanup
rm get-pip.py >& /dev/null
#rm -rf ~/.local &> /dev/null

# Exit with the status from the python testing step
exit ${err:?}

