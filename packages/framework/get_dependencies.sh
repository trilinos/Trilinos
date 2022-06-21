#!/bin/bash 
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# Data that needs to be updated when GenConfig changes!
genconfig_sha1=924a08af66f0a0573b5dd1128179731489339aec

# The following code contains no changing data

pushd $PWD &> /dev/null

cd ${script_dir}

function tril_genconfig_clone_or_update_repo() {
  git_url=$1
  sub_dir=$2
  has_submodules=$3
  head_sha=$4

  pushd $PWD &> /dev/null

  echo

  if [[ -d ${sub_dir} ]] ; then
    echo "STATUS: ${sub_dir}: Fetching remote repo"
    cd ${sub_dir}
    git fetch
  else
    echo "STATUS: ${sub_dir}: Cloning from '${git_url}'"
    git clone ${git_url} ${sub_dir}
    cd ${sub_dir}
  fi

  if [[ ! -z ${head_sha} ]]; then
    git checkout -f ${head_sha}
  else
    echo "STATUS: ${sub_dir}: Merging tip of remote tracking branch"
    git merge @{u}
  fi

  if [[ "${has_submodules}" == "has-submodules" ]] ; then
    echo
    echo "STATUS: ${sub_dir}: Update submodules"
    git submodule update --force --init --recursive
    cd - > /dev/null
  elif [[ "${has_submodules}" != "" ]] ; then
    echo "ERROR: argument '${has_submodules}' not allowed!  Only 'has-submodules' or ''!"
    exit 1
  fi

  popd &> /dev/null
}

# Clone or update the repos

tril_genconfig_clone_or_update_repo \
  git@gitlab-ex.sandia.gov:trilinos-devops-consolidation/code/GenConfig.git \
  GenConfig  has-submodules ${genconfig_sha1}

if [[ "$ini_file_option" == "--srn" ]] ; then
  tril_genconfig_clone_or_update_repo \
    git@cee-gitlab.sandia.gov:trilinos-project/srn-ini-files.git \
    srn-ini-files
elif [[ "$ini_file_option" == "--son" ]] ; then
  tril_genconfig_clone_or_update_repo \
    git@gitlab-ex.sandia.gov:trilinos-project/son-ini-files.git \
    son-ini-files
elif [[ "$ini_file_option" != "" ]] ; then
  echo "ERROR: Option '${ini_file_option}' not allowed! Must select '--son', '--srn' or ''."
  exit 1
fi

# Set up symlinks to the desired *.ini files
echo
cd ${script_dir}/GenConfig/deps/LoadEnv/ini_files
if [[ -d ${script_dir}/srn-ini-files ]] && [[ "$ini_file_option" == "--srn" ]]; then
    echo "STATUS: Link files from srn-ini-files"
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
elif [[ -d ${script_dir}/son-ini-files ]] && [[ "$ini_file_option" == "--son" ]]; then
    echo "STATUS: Link files from son-ini-files"
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${script_dir}/son-ini-files/trilinos/framework/supported-systems.ini
else
    echo "STATUS: Link files from ini-files"
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
echo
echo "INFO: You selected the following LoadEnv ini files:"
echo
find GenConfig/deps/LoadEnv/ini_files -type l -exec ls -lta {} \;
echo
echo "INFO: You selected the following GenConfig ini files:"
echo
find  GenConfig/ini_files -type l -exec ls -lta {} \;
echo
echo "INFO: If these symlinks do not point to the desired ini files, please re-run:"
echo
echo "    $0 [--son|--srn]"
popd > /dev/null
