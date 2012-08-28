INCLUDE(TribitsSubPackageMacros)
INCLUDE(TribitsLibraryMacros)


TRIBITS_SUBPACKAGE(commsplitter)

SET(HEADERS "")
SET(NONINSTALLHEADERS "")
SET(SOURCES "")

APPEND_SET(NOINSTALLHEADERS
  commsplitter.h
)

APPEND_SET(SOURCES
  init.c
  logging.c
  wrappers.c
)

############# libcommsplitter ##########################

TRIBITS_ADD_LIBRARY(
   commsplitter
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${SOURCES}
   DEPLIBS ${DEPLIBS} -lz
)

TRIBITS_SUBPACKAGE_POSTPROCESS()
