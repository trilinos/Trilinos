# This demo builds all sources into a library, and then compiles
# and links  main  against the library and all Trilinos libraries.

SET(SOURCES
	src_file.cpp 
  )

SET(HEADERS
	src_file.hpp
  )

INCLUDE_DIRECTORIES ( ./  ./src_subdirTBD
                     ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

ADD_LIBRARY(myappLib ${SOURCES} ${HEADERS})

ADD_EXECUTABLE(MyApp.exe main_file.cpp)

TARGET_LINK_LIBRARIES(MyApp.exe  myappLib  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}) 
