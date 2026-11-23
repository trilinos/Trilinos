TRIBITS_SUBPACKAGE(Exodus_cpp)

if (SEACASExodus_ENABLE_STATIC)
INCLUDE(InstallSymLink)
endif()

TRIBITS_INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

SET(SOURCES "")
SET(HEADERS "")

APPEND_GLOB(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.C)
APPEND_GLOB(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_LIBRARY(
  exodus_cpp
  HEADERS ${HEADERS}
  SOURCES ${SOURCES} ${DEP_SOURCES}
)

set_property(TARGET exodus_cpp PROPERTY C_STANDARD 99)
set_property(TARGET exodus_cpp PROPERTY C_EXTENSIONS ON)

if (SEACASExodus_ENABLE_SHARED)
   if (NOT BUILD_SHARED_LIBS)
      ADD_LIBRARY(
              exodus_cpp_shared
              SHARED
              ${SOURCES} ${DEP_SOURCES}
      )
      target_include_directories(exodus_cpp_shared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
      target_link_libraries(exodus_cpp_shared PUBLIC SEACASExodus::all_libs suplib_cpp)
      target_link_libraries(exodus_cpp_shared PUBLIC Netcdf::all_libs)
      set_property(TARGET exodus_cpp_shared PROPERTY C_STANDARD 99)
      # This keeps the library out of the `all_libs` targets...
      set_target_properties(exodus_cpp_shared PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
      set_target_properties(exodus_cpp_shared PROPERTIES OUTPUT_NAME exodus_cpp)
      INSTALL(TARGETS exodus_cpp_shared DESTINATION ${CMAKE_INSTALL_LIBDIR})
   endif()
endif()

# Many users of the SEACAS installations still want a static exodus library...
# But if only building a static library, it will be added in the `TRIBITS_ADD_LIBRARY` call.
if (SEACASExodus_ENABLE_STATIC AND SEACASExodus_ENABLE_SHARED)
   # NOTE: If use `TRIBITS_ADD_LIBRARY` here, then both exodus_cpp and exodus_cpp_static linked to all executables.
   ADD_LIBRARY(
           exodus_cpp_static
           STATIC
           ${SOURCES} ${DEP_SOURCES}
   )
   target_include_directories(exodus_cpp_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
   target_link_libraries(exodus_cpp_static PUBLIC SEACASExodus::all_libs suplib_cpp)
   target_link_libraries(exodus_cpp_static PUBLIC Netcdf::all_libs)
   set_property(TARGET exodus_cpp_static PROPERTY C_STANDARD 99)
   # This keeps the library out of the `all_libs` targets...
   set_target_properties(exodus_cpp_static PROPERTIES TRIBITS_TESTONLY_LIB TRUE)
   set_target_properties(exodus_cpp_static PROPERTIES OUTPUT_NAME exodus_cpp)
   INSTALL(TARGETS exodus_cpp_static DESTINATION ${CMAKE_INSTALL_LIBDIR})

   # For backwards compatibility, the C++ exodus library was called libexoIIv2cpp.a
   # Create a symbolic link from libexodus.a to libexoIIv2cpp.a -- do only for static lib
   InstallSymLink(libexodus_cpp.a ${CMAKE_INSTALL_PREFIX}/lib/libexoIIv2cpp.a)
endif()

TRIBITS_SUBPACKAGE_POSTPROCESS()
