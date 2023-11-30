FILE(GLOB PYROL_SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
MESSAGE("PYROL_SRC = ${PYROL_SRC}")
FILE(COPY ${PYROL_SRC} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

MESSAGE("CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")

file(STRINGS ${CMAKE_CURRENT_BINARY_DIR}/../binder/pyrol.sources BINDER_SRCS)
list(TRANSFORM BINDER_SRCS PREPEND "${CMAKE_CURRENT_BINARY_DIR}/../binder/")
list(APPEND PYROL_SRC ${BINDER_SRCS})

MESSAGE("PYROL_SRC with binder = ${PYROL_SRC}")

pybind11_add_module(pyrol ${PYROL_SRC})
target_include_directories(pyrol PUBLIC ${Mpi4Py_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/../binder "$<BUILD_INTERFACE:${ROL_all_include_dirs}>")
target_compile_features(pyrol PUBLIC cxx_std_14)

foreach(depPkg IN LISTS ROL_LIB_ENABLED_DEPENDENCIES) 
    target_link_libraries(pyrol PUBLIC ${depPkg}::all_libs)
endforeach()
target_link_libraries(pyrol PUBLIC ${Trilinos_EXTRA_LINK_FLAGS})
set_target_properties(pyrol PROPERTIES SUFFIX ".so")

INSTALL(TARGETS pyrol
  DESTINATION ${PyROL_INSTALL_DIR})

add_custom_command(TARGET pyrol POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/pyrol.so ${CMAKE_CURRENT_BINARY_DIR}/../pyrol/.
    COMMENT "Copy ${PROJECT_BINARY_DIR}/src/PyROL.so"
)
