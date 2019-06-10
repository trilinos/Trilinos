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
