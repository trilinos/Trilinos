#!/bin/bash 
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

pushd $PWD

cd ${script_dir}

# Update GenConfig submodules
# NOTE: Add --remote below to pull in the latest tip of the GenConfig tracked branch
git submodule update --init --recursive ${script_dir}/GenConfig

# Update ini file submodules
git submodule update --init --remote ${script_dir}/srn-ini-files || true
git submodule update --init --remote ${script_dir}/son-ini-files || true

# Point to srn ini files if they exist
if [[ ! -z "$(ls ./srn-ini-files)" && "$ini_file_option" == "--srn" ]]; then
    cd GenConfig/deps/LoadEnv/ini_files
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
    ln -sf ${script_dir}/ini-files/supported-envs.ini
elif [[ ! -z "$(ls ./son-ini-files)" && "$ini_file_option" == "--son" ]]; then
    cd GenConfig/deps/LoadEnv/ini_files
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/supported-systems.ini
    ln -sf ${script_dir}/ini-files/supported-envs.ini
else
    echo "ERROR: Unable to clone ini files repo or invalid ini_file_option=$ini_file_option"
    echo "Please ensure you have access to the ini files repo"
    echo ""
    echo "Usage: ./get_dependencies.sh [--srn,--son]"
    exit 1
fi

cd ../../../ini_files
ln -sf ${script_dir}/ini-files/config-specs.ini
ln -sf ${script_dir}/ini-files/supported-config-flags.ini

popd
