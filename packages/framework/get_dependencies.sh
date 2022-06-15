#!/bin/bash 
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

pushd $PWD

cd ${script_dir}

function tril_genconfig_clone_or_update_repo() {
  git_url=$1
  sub_dir=$2
  has_submodules=$3

  if [[ -e ${sub_dir} ]] ; then
    cd ${sub_dir}
    git pull
    cd - > /dev/null
  else
    git clone ${git_url} ${sub_dir}
  fi

  if [[ "${has_submodules}" == "has-submodules" ]] ; then
    cd ${sub_dir}
    git submodule update --init --recursive
    cd - > /dev/null
  elif [[ "${has_submodules}" != "" ]] ; then
    echo "ERROR: argument '${has_submodules}' not allowed!  Only 'has-submodules' or ''!"
    exit 1
  fi
}

# Clone or update the repos

tril_genconfig_clone_or_update_repo \
  git@gitlab-ex.sandia.gov:trilinos-devops-consolidation/code/GenConfig.git \
  GenConfig  has-submodules

if [[ "$ini_file_option" == "--srn" ]] ; then
  tril_genconfig_clone_or_update_repo \
    git@cee-gitlab.sandia.gov:trilinos-project/srn-ini-files.git \
    srn-ini-files
fi

if [[ "$ini_file_option" == "--son" ]] ; then
  tril_genconfig_clone_or_update_repo \
    git@gitlab-ex.sandia.gov:trilinos-project/son-ini-files.git \
    son-ini-files
fi

# Set up symlinks to the desired *.ini files
cd GenConfig/deps/LoadEnv/ini_files
if [[ -d ${script_dir}/srn-ini-files ]] && [[ "$ini_file_option" == "--srn" ]]; then
    echo "Link files from srn-ini-files"
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
    ln -sf ${script_dir}/ini-files/supported-envs.ini
elif [[ -d ${script_dir}/son-ini-files ]] && [[ "$ini_file_option" == "--son" ]]; then
    echo "Link files from son-ini-files"
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/supported-systems.ini
    ln -sf ${script_dir}/ini-files/supported-envs.ini
else
    echo "Link files from init-files"
    [ -e ${script_dir}/ini-files/environment-specs.ini ] && ln -sf ${script_dir}/ini-files/environment-specs.ini
    [ -e ${script_dir}/ini-files/supported-systems.ini ] && ln -sf ${script_dir}/ini-files/supported-systems.ini
    ln -sf ${script_dir}/ini-files/supported-envs.ini
fi
cd - > /dev/null
cd GenConfig/ini_files

ln -sf ${script_dir}/ini-files/config-specs.ini
ln -sf ${script_dir}/ini-files/supported-config-flags.ini

popd > /dev/null
