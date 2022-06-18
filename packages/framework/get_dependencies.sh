#!/bin/bash 
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

pushd $PWD

cd ${script_dir}

function tril_genconfig_clone_or_update_repo() {
  pushd $PWD
  git_url=$1
  sub_dir=$2
  has_submodules=$3
  head_sha=$4

  if [[ -e ${sub_dir} ]] ; then
    cd ${sub_dir}
    git fetch
  else
    git clone ${git_url} ${sub_dir}
    cd ${sub_dir}
  fi

  if [[ ! -z ${head_sha} ]]; then
    git checkout -f ${head_sha}
  fi

  if [[ "${has_submodules}" == "has-submodules" ]] ; then
    git submodule update --force --init --recursive
    cd - > /dev/null
  elif [[ "${has_submodules}" != "" ]] ; then
    echo "ERROR: argument '${has_submodules}' not allowed!  Only 'has-submodules' or ''!"
    exit 1
  fi
  popd > /dev/null
}

# Clone or update the repos

tril_genconfig_clone_or_update_repo \
  git@gitlab-ex.sandia.gov:trilinos-devops-consolidation/code/GenConfig.git \
  GenConfig  has-submodules b47cbdb1145a9d13f5bc8d4e7e691bab1dcdb60d

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
cd ${script_dir}/GenConfig/deps/LoadEnv/ini_files
if [[ -d ${script_dir}/srn-ini-files ]] && [[ "$ini_file_option" == "--srn" ]]; then
    echo "Link files from srn-ini-files"
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
elif [[ -d ${script_dir}/son-ini-files ]] && [[ "$ini_file_option" == "--son" ]]; then
    echo "Link files from son-ini-files"
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/supported-systems.ini
else
    echo "Link files from ini-files"
    [ -e ${script_dir}/ini-files/environment-specs.ini ] && ln -sf ${script_dir}/ini-files/environment-specs.ini
    [ -e ${script_dir}/ini-files/supported-systems.ini ] && ln -sf ${script_dir}/ini-files/supported-systems.ini
fi
ln -sf ${script_dir}/ini-files/supported-envs.ini
cd - > /dev/null

cd ${script_dir}/GenConfig/ini_files
ln -sf ${script_dir}/ini-files/config-specs.ini
ln -sf ${script_dir}/ini-files/supported-config-flags.ini

# Print summary of ini file settings
cd ${script_dir}
echo "You selected the following LoadEnv ini files:"
ls -lat GenConfig/deps/LoadEnv/ini_files
echo ""
echo "You selected the following GenConfig ini files:"
ls -lat GenConfig/ini_files
echo "If these symlinks do not point to the desired ini files, please re-run:"
echo "    $0 [--son, --srn]"
popd > /dev/null
