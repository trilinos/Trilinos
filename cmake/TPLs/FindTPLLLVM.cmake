find_package(LLVM CONFIG)

IF (NOT LLVM_FOUND)
  MESSAGE(FATAL_ERROR "LLVM was not found!")
ENDIF()

message(STATUS "Using LLVM includes in: ${LLVM_INCLUDE_DIRS}")
