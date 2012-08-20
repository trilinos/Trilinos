############# Function to generate XDR source and header Files ###############

function (TriosProcessXDR path)

   GET_FILENAME_COMPONENT(file ${path} NAME_WE)

   add_custom_command(
     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.c
     COMMAND rpcgen -Cc ${path}
        | sed -e "\"s#include.*${file}.*#include <${file}.h>#\""
        > ${CMAKE_CURRENT_BINARY_DIR}/${file}.c
     DEPENDS ${path} ${file}.h)


   add_custom_command(
     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.h
     COMMAND rpcgen -Ch ${path}
        | sed -e "s#rpc/rpc.h#${file}.h#"
        > ${CMAKE_CURRENT_BINARY_DIR}/${file}.h
     DEPENDS ${path})

   # Need target to force construction of nssi_types_xdr.{c,h} and nnti_xdr.{c,h}
   add_custom_target(${file} ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${file}.c ${CMAKE_CURRENT_BINARY_DIR}/${file}.h)

   INCLUDE_DIRECTORIES(BEFORE SYSTEM ${CMAKE_CURRENT_BINARY_DIR})

   # These Flags convert the "extra comma" error into a warning when using the -pedantic flag
   IF (CMAKE_COMPILER_IS_GNUCXX)
     #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
     SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
     SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -isystem ${CMAKE_CURRENT_BINARY_DIR}")
   ENDIF (CMAKE_COMPILER_IS_GNUCXX)

endfunction (TriosProcessXDR)

