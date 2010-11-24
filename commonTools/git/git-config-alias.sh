# Set up some shortcut commands

# Note that the dates are shown in the format:
#
#   [absolute author date] (relative committer date)
#
# The author date is when the commit was first created by the author.  The
# committer date is the date of the last time this commit was changed in
# anyway (e.g. in a rebase or an amend).  These dates are often the same but
# not if a rebase occurs.  Anyway, I think showing both of these dates is
# helpful.

# Shorter versions of log
eg config --global alias.log-short "log --pretty=format:'%Cgreen%h%Creset \"%s\"%nAuthor: %an <%ae>%nDate:   %ad (%cr)%n'"
eg config --global alias.log-oneline "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)'"
eg config --global alias.log-oneline-nc "log --pretty=format:'%h \"%s\" <%ae> [%ad] (%cr)'"
eg config --global alias.log-local "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)' origin.."

# Summarizing changes locally and remotely
eg config --global alias.local-stat "!git status ; echo ; echo 'Commits in local repo not yet pushed to origin/master:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)' --name-status HEAD ^origin/master"
eg config --global alias.remote-stat "!git status ; echo ; echo 'Commits in origin/master not in local repo:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)' --name-status ^HEAD origin/master"
eg config --global alias.local-stat-short "!git status ; echo ; echo 'Commits in local repo not yet pushed to origin/master:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)' --shortstat --dirstat=0 HEAD ^origin/master"
eg config --global alias.remote-stat-short "!git status ; echo ; echo 'Commits in origin/master not in local repo:' ; echo ; git log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%cr)' --shortstat --dirstat=0 ^HEAD origin/master"

# Diffing text files like latex (undocumented option --color-words)
eg config --global alias.wlog "log --color-words"
eg config --global alias.wdiff "diff --color-words"
