# Set up some shortcut commands
eg config --global alias.log-short "log --pretty=format:'%Cgreen%h%Creset \"%s\"%nAuthor: %an <%ae>%nDate:   %ad (%ar)%n'"
eg config --global alias.log-oneline "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)'"
eg config --global alias.log-oneline-nc "log --pretty=format:'%h \"%s\" <%ae> [%ad] (%ar)'"
eg config --global alias.log-local "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' origin.."
eg config --global alias.local-stat "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' --name-status HEAD ^origin/master"
eg config --global alias.remote-stat "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> [%ad] (%ar)' --name-status ^HEAD origin/master"
