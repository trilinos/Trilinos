INCLUDE(MessageWrapper)
INCLUDE(Split)

#
# @FUNCTION: TRIBITS_DIR_IS_BASEDIR()
#
# Function to determine if a given path is a base dir of another path.
#
# Usage::
#
#   TRIBITS_DIR_IS_BASEDIR(<absBaseDir> <absFullDir> <isBaseDirVarOut>)
#
# If the absolute path ``<absBaseDir>`` is a subdir of the absolute path
# ``<absFullDir>``, then the variable ``<isBaseDirVarOut>`` is set to
# ``TRUE``.  Otherwise, ``<isBaseDirVarOut>`` is set to ``FALSE``.
#
# For example, the output var ``isBaseDir`` would be set to ``TRUE`` in the
# following examples::
#
#   TRIBITS_DIR_IS_BASEDIR(/some/base/path /some/base/path/more isBaseDir)
#
#   TRIBITS_DIR_IS_BASEDIR(/some/base/path /some/base/path isBaseDir)
#
# However, in the following examples, ``isBaseDir`` would be set to ``FALSE``::
#
#   TRIBITS_DIR_IS_BASEDIR(/some/base/path/more /some/base/path isBaseDir)
#
#   TRIBITS_DIR_IS_BASEDIR(/some/base/path /some/other/path isBaseDir)
#
FUNCTION(TRIBITS_DIR_IS_BASEDIR  absBaseDir  absFullDir  isBaseDirVarOut)

  # Assume not base dir by default unless we find it is
  SET(isBaseDir FALSE)

  STRING(LENGTH "${absBaseDir}" absBaseDirLen)
  STRING(LENGTH "${absFullDir}" absFullDirLen)

  IF (absBaseDir STREQUAL absFullDir)
    SET(isBaseDir TRUE)
  ELSEIF (NOT absBaseDirLen GREATER absFullDirLen)
    STRING(FIND "${absFullDir}" "${absBaseDir}/" baseDirIdx)
    IF (baseDirIdx EQUAL 0)
      SET(isBaseDir TRUE)
    ENDIF()
  ENDIF()

  SET(${isBaseDirVarOut} ${isBaseDir} PARENT_SCOPE)

ENDFUNCTION()


#
# @FUNCTION: TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR()
#
# Returns the array of directories below a base directory.
#
# Usage::
#
#  TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR(<absBaseDir> <absFullDir>
#     <trailingDirArrayVarOut>)
#
# The following show examples of what this returns:
#
#   TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR("/a/b/c" "/a/b/c", dirArray)
#     => dirArray = ""
#
#   TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR("/a/b/c" "/a/b/c/d", dirArray)
#     => dirArray = "d"
#
#   TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR("/a/b/c" "/a/b/c/d/e", dirArray)
#     => dirArray = "d;e"
#
FUNCTION(TRIBITS_GET_DIR_ARRAY_BELOW_BASE_DIR  absBaseDir  absFullDir
  trailingDirArrayVarOut
  )

  TRIBITS_DIR_IS_BASEDIR("${absBaseDir}" "${absFullDir}" isBaseDir)
  IF (NOT isBaseDir)
    MESSAGE_WRAPPER(FATAL_ERROR
      "ERROR: '${absBaseDir}' is not a base dir of '${absFullDir}'")
  ENDIF()

  STRING(LENGTH "${absBaseDir}" absBaseDirLen)
  STRING(LENGTH "${absFullDir}" absFullDirLen)

  IF (absBaseDirLen EQUAL absFullDirLen)
    SET(trailingDirArray "")
  ELSE()
    MATH(EXPR trailingDirsStrStartIdx "${absBaseDirLen}+1")
    STRING(SUBSTRING "${absFullDir}" ${trailingDirsStrStartIdx} -1 trailingDirsStr)
    SPLIT("${trailingDirsStr}" "/" trailingDirArray)
  ENDIF()

  SET(${trailingDirArrayVarOut} "${trailingDirArray}" PARENT_SCOPE)

ENDFUNCTION()
