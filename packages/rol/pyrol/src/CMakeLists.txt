FILE(GLOB PYROL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

# The generated sources: pyrol.cpp alone with PYROL_BINDER_USE_ONE_FILE, and
# otherwise the files Binder lists in pyrol.sources.
SET(PYROL_BINDER_DIR ${CMAKE_CURRENT_BINARY_DIR}/../binder)
IF(PYROL_USE_ONE_FILE)
  SET(PYROL_BINDER_SRC pyrol.cpp)
ELSE()
  FILE(STRINGS ${PYROL_BINDER_DIR}/pyrol.sources PYROL_BINDER_SRC)
ENDIF()
list(TRANSFORM PYROL_BINDER_SRC PREPEND "${PYROL_BINDER_DIR}/")
list(APPEND PYROL_SRC ${PYROL_BINDER_SRC})

pybind11_add_module(pyrol ${PYROL_SRC})
target_include_directories(pyrol PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PYROL_BINDER_DIR} "$<BUILD_INTERFACE:${ROL_all_include_dirs}>")
target_compile_features(pyrol PUBLIC cxx_std_17)

foreach(depPkg IN LISTS ROL_LIB_ENABLED_DEPENDENCIES)
    target_link_libraries(pyrol PUBLIC ${depPkg}::all_libs)
endforeach()
target_link_libraries(pyrol PUBLIC ${Trilinos_EXTRA_LINK_FLAGS})
set_target_properties(pyrol PROPERTIES SUFFIX ".so")

INSTALL(TARGETS pyrol
  DESTINATION ${PyROL_INSTALL_DIR})

add_custom_command(TARGET pyrol POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/pyrol.so ${CMAKE_CURRENT_BINARY_DIR}/../pyrol/.
    COMMENT "Copy pyrol.so into the build tree's pyrol package"
)
