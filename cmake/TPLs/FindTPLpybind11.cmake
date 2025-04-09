find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11)

IF (NOT pybind11_FOUND)
  MESSAGE(FATAL_ERROR "pybind11 was not found!")
ENDIF()
