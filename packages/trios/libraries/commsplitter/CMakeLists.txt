
INCLUDE(PackageLibraryMacros)

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

PACKAGE_ADD_LIBRARY(
   commsplitter
   HEADERS ${HEADERS}
   NOINSTALLHEADERS ${NOINSTALLHEADERS}
   SOURCES ${SOURCES}
   DEPLIBS ${DEPLIBS}
)

