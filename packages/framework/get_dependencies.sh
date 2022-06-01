#!/bin/bash
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

pushd $PWD

cd ${script_dir}

# Update submodules
git submodule update --init --recursive --remote || true

# Point to srn ini files if they exist
if [[ ! -z "$(ls ./srn-ini-files)" && "$ini_file_option" == "--srn" ]]; then
    cd GenConfig/deps/LoadEnv/ini_files
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
else
    cd GenConfig/deps/LoadEnv/
    git reset --hard HEAD
fi

popd
