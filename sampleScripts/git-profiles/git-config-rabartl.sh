# Set your username
eg config --global user.name "Roscoe A. Bartlett"
# Set your email address
eg config --global user.email rabartl@sandia.gov
# Use colorized output when it makes sense
eg config --global color.ui true
# Set up some shortcut commands
eg config --global alias.log-short1 "log --pretty=format:'%Cgreen%H%Creset%nAuthor: %an <%ae>%nDate: %cd (%ar)%n%n%s%n'"
eg config --global alias.log-short2 "log --pretty=format:'%Cgreen%h%Creset \"%s\"%nAuthor: %an <%ae>%nDate: %ad (%ar)%n'"
eg config --global alias.log-online1 "log --pretty=format:'%Cgreen%h%Creset \"%s\" <%ae> %ad'"
