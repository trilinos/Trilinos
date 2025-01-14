#!/bin/bash -e
ini_file_option=$1
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# Data that needs to be updated when GenConfig changes!
genconfig_sha1=88c44e347c0377a170ec9ca45a47732a9630b4ec


# The following code contains no changing data

pushd $PWD &> /dev/null

cd ${script_dir}

function tril_genconfig_assert_pwd_is_git_repo() {
  if [[ ! -d "${PWD}/.git" ]] ; then
    echo "ERROR: The directrory ${PWD} is not a git repository!"
    echo "  ==> Please delete ${PWD} and run again"
    exit 1
  fi
}

function retry_command() {
  cmd=$1
  ${cmd} || { echo "Retrying after 1m..." ; sleep 60 ; ${cmd} ; } || { echo "Retrying after 5m..." ; sleep 300 ; ${cmd} ; }
}

function tril_genconfig_clone_or_update_repo() {
  git_url=$1
  sub_dir=$2
  has_submodules=$3
  head_sha=$4

  pushd $PWD &> /dev/null

  echo

  if [[ -d ${sub_dir} ]] ; then
    cd ${sub_dir}
    # Validate correct remote and abort if not correct
    remote=$(git remote get-url origin)
    normalized_remote=$(python3 ${script_dir}/normalize_git_repo_url.py ${remote})
    normalized_git_url=$(python3 ${script_dir}/normalize_git_repo_url.py ${git_url})
    if [[ "${normalized_remote}" != "${normalized_git_url}" ]] ; then
      echo "ERROR: Current remote origin does not match expected!" >&2
      echo "   Current Remote:  ${remote}" >&2
      echo "   Expected Remote: ${git_url}" >&2
      echo "" >&2
      echo "Please remove/move '$(pwd)' and re-run this script" >&2
      echo "" >&2
      echo "" >&2
      exit 1
    fi
    # Update remote repo (which points to correct remote)
    echo "STATUS: ${sub_dir}: Fetching remote repo"
    tril_genconfig_assert_pwd_is_git_repo
    cmd="git fetch"
    retry_command "${cmd}"
  else
    echo "STATUS: ${sub_dir}: Cloning from '${git_url}'"
    cmd="git clone ${git_url} ${sub_dir}"
    retry_command "${cmd}"
    cd ${sub_dir}
  fi

  if [[ ! -z ${head_sha} ]]; then
    tril_genconfig_assert_pwd_is_git_repo
    git checkout -f ${head_sha}
  else
    echo "STATUS: ${sub_dir}: Merging tip of remote tracking branch"
    tril_genconfig_assert_pwd_is_git_repo
    git merge @{u}
  fi

  if [[ "${has_submodules}" == "has-submodules" ]] ; then
    echo
    echo "STATUS: ${sub_dir}: Update submodules"
    cmd="git submodule update --force --init --recursive"
    retry_command "${cmd}"
    cd - > /dev/null
  elif [[ "${has_submodules}" != "" ]] ; then
    echo "ERROR: argument '${has_submodules}' not allowed!  Only 'has-submodules' or ''!"
    exit 1
  fi

  popd &> /dev/null
}

# Clone GenConfig from GitHub
tril_genconfig_clone_or_update_repo \
  https://github.com/sandialabs/GenConfig.git \
  GenConfig  has-submodules ${genconfig_sha1}

if [[ "$ini_file_option" == "--srn" ]] ; then
  #Clone srn-ini-files from cee-gitlab
  tril_genconfig_clone_or_update_repo \
    git@cee-gitlab.sandia.gov:trilinos-project/srn-ini-files.git \
    srn-ini-files
  
elif [[ "$ini_file_option" == "--son" ]] ; then
  #Clone son-ini-files from gitlab-ex
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
relative_path_to_script_dir="../../../.."
if [[ -d ${script_dir}/srn-ini-files ]] && [[ "$ini_file_option" == "--srn" ]]; then
    echo "STATUS: Link files from srn-ini-files"
    ln -sf ${relative_path_to_script_dir}/srn-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${relative_path_to_script_dir}/srn-ini-files/trilinos/framework/supported-systems.ini
elif [[ -d ${script_dir}/son-ini-files ]] && [[ "$ini_file_option" == "--son" ]]; then
    echo "STATUS: Link files from son-ini-files"
    ln -sf ${relative_path_to_script_dir}/son-ini-files/trilinos/framework/environment-specs.ini
    ln -sf ${relative_path_to_script_dir}/son-ini-files/trilinos/framework/supported-systems.ini
else
    echo "STATUS: Link files from ini-files"
    [ -e ${relative_path_to_script_dir}/ini-files/environment-specs.ini ] && ln -sf ${relative_path_to_script_dir}/ini-files/environment-specs.ini
    [ -e ${relative_path_to_script_dir}/ini-files/supported-systems.ini ] && ln -sf ${relative_path_to_script_dir}/ini-files/supported-systems.ini
fi
ln -sf ${relative_path_to_script_dir}/ini-files/supported-envs.ini
cd - > /dev/null

cd ${script_dir}/GenConfig/ini_files
relative_path_to_script_dir="../.."
ln -sf ${relative_path_to_script_dir}/ini-files/config-specs.ini
ln -sf ${relative_path_to_script_dir}/ini-files/supported-config-flags.ini

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
