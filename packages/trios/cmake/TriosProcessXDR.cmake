############# Function to generate XDR source and header Files ###############

function (TriosProcessXDR path)

   GET_FILENAME_COMPONENT(file ${path} NAME_WE)

   add_custom_command(
     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.c
     COMMAND rpcgen -Cc ${path}
        | sed -e "\"s#include.*${file}.*#include <${file}.h>#\""
        | grep -v \"register int32_t \\*buf\;\"
        > ${CMAKE_CURRENT_BINARY_DIR}/${file}.c
     DEPENDS ${path} ${file}.h)


   add_custom_command(
     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${file}.h
     COMMAND rpcgen -Ch ${path}
        | sed -e "\"s#rpc/rpc.h#${file}.h#\""
        | perl -pe \"BEGIN{undef $$/\;} s/\(enum\\s\\w+\\s\\{\\n\(\\s*.*?,\\n\)*?\\s*.*?\),\(\\n\\s*\\}\;\)/\\1\\3/smg\"
        > ${CMAKE_CURRENT_BINARY_DIR}/${file}.h
     DEPENDS ${path})

   # Need target to force construction of nssi_types_xdr.{c,h} and nnti_xdr.{c,h}
   add_custom_target(${file} ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${file}.c ${CMAKE_CURRENT_BINARY_DIR}/${file}.h)

endfunction (TriosProcessXDR)
