find_package(Git)

if (Git_FOUND)

  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --dirty --always --match=""
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE TRILINOS_GIT_SHA
    RESULT_VARIABLE RESULT
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(RESULT EQUAL 0)
    message(STATUS "Current Git SHA: ${TRILINOS_GIT_SHA}")
  else()
    message(STATUS "Git returned code \"${RESULT}\", message \"${TRILINOS_GIT_SHA}\"")
    set(TRILINOS_GIT_SHA "UNDEFINED")
  endif()
else()
  message(STATUS "Git not found. Could not collect SHA.")
  set(TRILINOS_GIT_SHA "UNDEFINED")
endif()
configure_file(${Trilinos_SOURCE_DIR}/cmake/Trilinos_git_sha.h.in Trilinos_git_sha.h)
