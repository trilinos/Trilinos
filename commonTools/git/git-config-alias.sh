# Set up some shortcut commands

# Shorter versions of log
eg config --global alias.log-short "log --pretty=format:'%Cgreen%h%Creset \"%s\"%nAuthor: %an <%ae>%nDate:   %ad (%ar)%n'"
eg config --global alias.log-oneline "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)'"
eg config --global alias.log-oneline-nc "log --pretty=format:'%h \"%s\" <%ae> [%ad] (%ar)'"
eg config --global alias.log-local "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' origin.."

# Summarizing changes locally and remotely
eg config --global alias.local-stat "!git status ; echo ; echo 'Commits in local repo not yet pushed to origin/master:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' --name-status HEAD ^origin/master"
eg config --global alias.remote-stat "!git status ; echo ; echo 'Commits in origin/master not in local repo:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' --name-status ^HEAD origin/master"

# Diffing text files like latex (undocumented option --color-words)
eg config --global alias.wlog "log --color-words"
eg config --global alias.wdiff "diff --color-words"
