include(MessageWrapper)
include(Split)


# @FUNCTION: tribits_dir_is_basedir()
#
# Function to determine if a given path is a base dir of another path.
#
# Usage::
#
#   tribits_dir_is_basedir(<absBaseDir> <absFullDir> <isBaseDirVarOut>)
#
# If the absolute path ``<absBaseDir>`` is a subdir of the absolute path
# ``<absFullDir>``, then the variable ``<isBaseDirVarOut>`` is set to
# ``TRUE``.  Otherwise, ``<isBaseDirVarOut>`` is set to ``FALSE``.
#
# For example, the output var ``isBaseDir`` would be set to ``TRUE`` in the
# following examples::
#
#   tribits_dir_is_basedir(/some/base/path /some/base/path/more isBaseDir)
#
#   tribits_dir_is_basedir(/some/base/path /some/base/path isBaseDir)
#
# However, in the following examples, ``isBaseDir`` would be set to ``FALSE``::
#
#   tribits_dir_is_basedir(/some/base/path/more /some/base/path isBaseDir)
#
#   tribits_dir_is_basedir(/some/base/path /some/other/path isBaseDir)
#
function(tribits_dir_is_basedir  absBaseDir  absFullDir  isBaseDirVarOut)

  # Assume not base dir by default unless we find it is
  set(isBaseDir FALSE)

  string(LENGTH "${absBaseDir}" absBaseDirLen)
  string(LENGTH "${absFullDir}" absFullDirLen)

  if (absBaseDir STREQUAL absFullDir)
    set(isBaseDir TRUE)
  elseif (NOT absBaseDirLen GREATER absFullDirLen)
    string(FIND "${absFullDir}" "${absBaseDir}/" baseDirIdx)
    if (baseDirIdx EQUAL 0)
      set(isBaseDir TRUE)
    endif()
  endif()

  set(${isBaseDirVarOut} ${isBaseDir} PARENT_SCOPE)

endfunction()


# @FUNCTION: tribits_get_dir_array_below_base_dir()
#
# Returns the array of directories below a base directory.
#
# Usage::
#
#  tribits_get_dir_array_below_base_dir(<absBaseDir> <absFullDir>
#     <trailingDirArrayVarOut>)
#
# The following show examples of what this returns:
#
#   tribits_get_dir_array_below_base_dir("/a/b/c" "/a/b/c", dirArray)
#     => dirArray = ""
#
#   tribits_get_dir_array_below_base_dir("/a/b/c" "/a/b/c/d", dirArray)
#     => dirArray = "d"
#
#   tribits_get_dir_array_below_base_dir("/a/b/c" "/a/b/c/d/e", dirArray)
#     => dirArray = "d;e"
#
function(tribits_get_dir_array_below_base_dir  absBaseDir  absFullDir
  trailingDirArrayVarOut
  )

  tribits_dir_is_basedir("${absBaseDir}" "${absFullDir}" isBaseDir)
  if (NOT isBaseDir)
    message_wrapper(FATAL_ERROR
      "ERROR: '${absBaseDir}' is not a base dir of '${absFullDir}'")
  endif()

  string(LENGTH "${absBaseDir}" absBaseDirLen)
  string(LENGTH "${absFullDir}" absFullDirLen)

  if (absBaseDirLen EQUAL absFullDirLen)
    set(trailingDirArray "")
  else()
    math(EXPR trailingDirsStrStartIdx "${absBaseDirLen}+1")
    string(SUBSTRING "${absFullDir}" ${trailingDirsStrStartIdx} -1 trailingDirsStr)
    split("${trailingDirsStr}" "/" trailingDirArray)
  endif()

  set(${trailingDirArrayVarOut} "${trailingDirArray}" PARENT_SCOPE)

endfunction()
