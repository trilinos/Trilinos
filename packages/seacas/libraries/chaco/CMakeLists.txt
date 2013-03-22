INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)

TRIBITS_SUBPACKAGE(Chaco)


SET(CHACO_VERSION_MAJOR "1")
SET(CHACO_VERSION_MINOR "0")
SET(CHACO_VERSION_PATCH "0")
SET(CHACO_VERSION "${CHACO_VERSION_MAJOR}.${CHACO_VERSION_MINOR}")
SET(CHACO_VERSION_FULL "${CHACO_VERSION}.${CHACO_VERSION_PATCH}")

INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/main"
  "${CMAKE_CURRENT_SOURCE_DIR}/util"
)

SET(HEADERS
  main/chaco.h
)

SET( chaco_SUBDIRS
     assign bpmatch coarsen connect eigen graph inertial
     input internal klspiff klvspiff main misc optimize
     refine_map refine_part submain symmlq tinvit util)

#- Add sources for a target
#
#  ADD_SOURCES(<target> <source1> [<source2> ...])
#
function(add_sources target)
  # define the <target>_SRCS properties if necessary
  get_property(prop_defined GLOBAL PROPERTY ${target}_SRCS DEFINED)
  if(NOT prop_defined)
    define_property(GLOBAL PROPERTY ${target}_SRCS
      BRIEF_DOCS "Sources for the ${target} target"
      FULL_DOCS "List of source files for the ${target} target")
  endif()
  # create list of sources (absolute paths)
  set(SRCS)
  foreach(src IN LISTS ARGN)
    if(NOT IS_ABSOLUTE "${src}")
      get_filename_component(src "${src}" ABSOLUTE)
    endif()
    list(APPEND SRCS "${src}")
  endforeach()
  # append to global property
  set_property(GLOBAL APPEND PROPERTY "${target}_SRCS" "${SRCS}")
endfunction()

# descend into sub-directories
FOREACH ( dir ${chaco_SUBDIRS} )
  add_subdirectory(${dir})
ENDFOREACH ( dir )


get_property(SOURCES GLOBAL PROPERTY chaco_SRCS)
TRIBITS_ADD_LIBRARY(
	chaco
	HEADERS	${HEADERS}
	SOURCES ${SOURCES}
	DEPLIBS m
)

TRIBITS_SUBPACKAGE_POSTPROCESS()
