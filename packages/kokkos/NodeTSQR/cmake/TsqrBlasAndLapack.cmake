foreach (_LIBNAME "BLAS" "LAPACK")
  set (_CURRENT_LIB "TPL_${_LIBNAME}_LIBRARIES")
  set (_CURRENT_LIB_VALUE "${${_CURRENT_LIB}}")
  if ("${_CURRENT_LIB_VALUE}" STREQUAL "")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (STATUS "${_CURRENT_LIB} not defined: building outside Trilinos")
    endif()
    set ("HAVE_${_CURRENT_LIB}" False)
  else ("${_CURRENT_LIB_VALUE}" STREQUAL "")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (STATUS "${_CURRENT_LIB} = ${_CURRENT_LIB_VALUE}: building inside Trilinos")
    endif()
    set ("HAVE_${_CURRENT_LIB}" True)
  endif ("${_CURRENT_LIB_VALUE}" STREQUAL "")
endforeach (_LIBNAME "BLAS" "LAPACK")

# Some BLAS and LAPACK implementations are unified into a single
# library.  Thus, we can only be sure that no LAPACK library was
# provided, if neither of these variables were set.
if (NOT HAVE_TPL_BLAS_LIBRARIES AND NOT HAVE_TPL_LAPACK_LIBRARIES)
#if (NOT ("${HAVE_TPL_BLAS_LIBRARIES}" STREQUAL "True") AND NOT ("${HAVE_TPL_LAPACK_LIBRARIES}" STREQUAL "True"))

  if ("${PKG_DIR}" STREQUAL "")
    set (PKG_DIR "$ENV{HOME}/pkg")
  endif ("${PKG_DIR}" STREQUAL "")
  if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
    message (STATUS "Setting TPL_{BLAS,LAPACK}_LIBRARIES here, since they were not set already")
    message (STATUS "Assuming libraries live in subdirectories of PKG_DIR = ${PKG_DIR}")
  endif()

  #set (BLAS_AND_LAPACK_SETUP mkl)
  #set (BLAS_AND_LAPACK_SETUP goto)
  set (BLAS_AND_LAPACK_SETUP default)

  # For Goto's BLAS, you should set the GOTO_NUM_THREADS environment
  # variable to 1 at runtime -- otherwise, bad thread-related things
  # will happen.  Goto's implementation manages a thread pool for
  # internal parallelism, but this can interact unpleasantly with
  # multithreading at higher levels of your application.  Similarly, for
  # Intel MKL, you should set BLAS_NUM_THREADS (or is it
  # OMP_NUM_THREADS?) to 1, to avoid similar problems.

  if ("${BLAS_AND_LAPACK_SETUP}" STREQUAL "goto")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (STATUS "-- Goto BLAS and custom LAPACK")
    endif()

    set (BLAS_LIB_DIR "${PKG_DIR}/GotoBLAS2-1.11p1")
    set (BLAS_LIBS goto2)

    if ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")
      set (FORTRAN_RUNTIME_LIBRARY gfortran)
    endif ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")

    set (LAPACK_LIB_DIR "${PKG_DIR}/lapack-3.2.1")

    # LAPACK depends on the Fortran runtime library corresponding to
    # the Fortran compiler used to build it.  If necessary,
    # FORTRAN_RUNTIME_LIBRARY should be defined via a -D argument to
    # CMake (preferably in a "do-configure" - type script).  If not
    # specified, it defaults to "".
    #
    # CMake doesn't like an empty space in a list of libraries, so we
    # define LAPACK_LIBS based on whether FORTRAN_RUNTIME_LIBRARY is
    # empty.
    if ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")
      set (LAPACK_LIBS "-llapack -ltmglib")
    else ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")
      set (LAPACK_LIBS "-llapack -ltmglib -l${FORTRAN_RUNTIME_LIBRARY}")
    endif ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")

    set (TPL_BLAS_LIBRARIES "-L${BLAS_LIB_DIR} ${BLAS_LIBS}")
    set (TPL_LAPACK_LIBRARIES "-L${LAPACK_LIB_DIR} ${LAPACK_LIBS}")
        
  elseif ("${BLAS_AND_LAPACK_SETUP}" STREQUAL "mkl")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (STATUS "-- Intel MKL BLAS and LAPACK")
    endif()

    set (MKL_LIB_PATH "/usr/mill/pkg/intel_mkl/lib/em64t/")
    set (BLAS_LIB_DIR ${MKL_LIB_PATH})
    set (LAPACK_LIB_DIR ${MKL_LIB_PATH})
    set (BLAS_LIBS "-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread")
    # LAPACK is mixed into the BLAS libraries.
    set (LAPACK_LIBS "")

    set (TPL_BLAS_LIBRARIES "-L${BLAS_LIB_DIR} ${BLAS_LIBS}")
    # This line differs slightly, because LAPACK_LIBS is empty, 
    # and CMake doesn't like an empty space here.
    set (TPL_LAPACK_LIBRARIES "-L${LAPACK_LIB_DIR}")

  elseif ("${BLAS_AND_LAPACK_SETUP}" STREQUAL "default")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (STATUS "-- Default system BLAS and custom LAPACK")
    endif()

    # Default setup is to use the system BLAS library with a
    # custom-built LAPACK library.
    set (BLAS_LIB_DIR "/usr/lib64")
    set (BLAS_LIBS "-lblas")

    # LAPACK depends on the Fortran runtime library corresponding to
    # the Fortran compiler used to build it.  FORTRAN_RUNTIME_LIBRARY
    # should be set in the "do-configure" script, or via a -D argument
    # when invoking CMake directly from the command line.

    # CMake doesn't like an empty space in a list of libraries, so we
    # define LAPACK_LIBS based on whether FORTRAN_RUNTIME_LIBRARY is
    # empty.
    set (LAPACK_LIB_DIR "${PKG_DIR}/lapack-3.2.1")
    if ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")
      set (LAPACK_LIBS "-llapack -ltmglib")
    else ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")
      set (LAPACK_LIBS "-llapack -ltmglib -l${FORTRAN_RUNTIME_LIBRARY}")
    endif ("${FORTRAN_RUNTIME_LIBRARY}" STREQUAL "")

    set (TPL_BLAS_LIBRARIES "-L${BLAS_LIB_DIR} ${BLAS_LIBS}")
    set (TPL_LAPACK_LIBRARIES "-L${LAPACK_LIB_DIR} ${LAPACK_LIBS}")

  else ("${BLAS_AND_LAPACK_SETUP}" STREQUAL "goto")
    if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
      message (FATAL_ERROR "No BLAS / LAPACK libraries specified")
    endif()

  endif ("${BLAS_AND_LAPACK_SETUP}" STREQUAL "goto")

else ()
  if (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
    message (STATUS "Using previously set TPL_{BLAS,LAPACK}_LIBRARIES")
  endif()
endif ()

IF (NOT ${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE)
  message (STATUS "* TPL_BLAS_LIBRARIES = ${TPL_BLAS_LIBRARIES}")
  message (STATUS "* TPL_LAPACK_LIBRARIES = ${TPL_LAPACK_LIBRARIES}")
ENDIF()


  
