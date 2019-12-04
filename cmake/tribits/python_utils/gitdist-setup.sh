# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Get the base dir for the sourced script
SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "SCRIPT_DIR = '$SCRIPT_DIR'"

existing_gitdist=`which gitdist 2> /dev/null`
if [[ "${existing_gitdist}" == "" ]] ; then
  echo "Setting alias gitdist=${SCRIPT_DIR}/gitdist"
  alias gitdist=${SCRIPT_DIR}/gitdist
fi

# Source this with bash to load useful env for using gitdist
alias gitdist-status="gitdist dist-repo-status"
alias gitdist-mod="gitdist --dist-mod-only"
alias gitdist-mod-status="gitdist --dist-mod-only dist-repo-status"
function gitdist-repo-versions {
  gitdist "$@" --dist-no-color log -1 --pretty=format:"%h [%ad] <%ae>%n%s" | grep -v "^$"
}

# Setup for completions for git command and gitdist options commands!
complete -o default -o nospace -F _git -W "dist-repo-status --dist-help --dist-use-git --dist-repos --dist-not-repos --dist-version-file --dist-version-file2 --dist-no-color --dist-debug --dist-no-opt --dist-mod-only --dist-legend" gitdist gitdist-mod
complete -o default -o nospace -W "--dist-use-git --dist-repos --dist-not-repos --dist-no-color --dist-debug --dist-no-opt --dist-mod-only" gitdist-repo-versions
