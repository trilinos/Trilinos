#
# Bash utility functions for generation of documentation.
#


#
# Copy a file from temp if it has changed.
# 
# Usage:
#
#   $ update_if_different  <base_file_name>  <tmpext>
#
# This function will copy from:
#
#   <base_file_name>.<tmpext>
#
# to:
#
#    <base_file_name>
#
# if <base_file_name>.<tmpext> diffs with <base_file_name>.  This is used
# avoid triggering a change in the timestamp of the file <base_file_name>.
#
# This isused to create *.tmp files and then we only copy to the files
# actually mentioned in the Makefile if they are differnet.  That way if the
# files don't change then the function will not copy over them and update
# their time stamp and the makefile will not rebuild the documentation.  We
# don't want to rebuild documentation if nothing changed.
#
function update_if_different {
  base_file_name=$1
  tmpext=$2
  if [ -e $base_file_name ] ; then
    changes=`diff $base_file_name.$tmpext $base_file_name`
    if [ "$changes" != "" ] ; then
       echo "Copy updated file $base_file_name" 
       cp $base_file_name.$tmpext $base_file_name
    fi
  else
    echo "Copy to non-existing file $base_file_name" 
    cp $base_file_name.$tmpext $base_file_name 
  fi
}


function generate_git_version_file {

  _TRIBITS_TAG_PREFIX=`cat ../../../tribits_tag_prefix.txt`

  if [ -e ../../../.git ] ; then
    echo
    echo "Generating git version"
    echo
    echo `git describe --match="$_TRIBITS_TAG_PREFIX*"` > TribitsGitVersion.txt
  else
    echo "$_TRIBITS_TAG_PREFIX.{Unknown version}" > TribitsGitVersion.txt
  fi

}
