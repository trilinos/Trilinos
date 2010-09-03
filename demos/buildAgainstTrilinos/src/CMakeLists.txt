# This demo builds all sources into a library, and then compiles
# links main against the library.

SET(SOURCES
	src_file.cpp 
  )

SET(HEADERS
	src_file.hpp
  )

INCLUDE_DIRECTORIES (${MYAPP_TRILINOS_DIR}/include 
                     ./  ./src_subdirTBD
                     ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES (${MYAPP_TRILINOS_DIR}/lib ${Trilinos_TPL_LIBRARY_DIRS})

ADD_LIBRARY(myappLib ${SOURCES} ${HEADERS})

ADD_EXECUTABLE(MyApp.exe main_file.cpp)

TARGET_LINK_LIBRARIES(MyApp.exe  myappLib  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 
