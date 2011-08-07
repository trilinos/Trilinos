# Set your username
eg config --global user.name "Roscoe A. Bartlett"
# Set your email address
eg config --global user.email bartlettra@ornl.gov
# Use colorized output when it makes sense
eg config --global color.ui true
# Set up some shortcut commands
_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
$_SCRIPT_DIR/../../commonTools/git/git-config-alias.sh
