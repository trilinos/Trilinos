# Remove a directory tree recursively if it exists, otherwise do nothing
#
# Usage:
#
#   cmake -D DIR_TO_REMOVE=<dir> -P remvoe_dir_if_exists.cmake
#

if (EXISTS "${DIR_TO_REMOVE}")
  message("Removing '${DIR_TO_REMOVE}' ...")
  file(REMOVE_RECURSE "${DIR_TO_REMOVE}") 
else()
  message("Dir '${DIR_TO_REMOVE}' does not exist!")
endif()
