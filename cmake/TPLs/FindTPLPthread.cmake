
set(USE_THREADS FALSE)

if (NOT TPL_Pthread_INCLUDE_DIRS AND NOT TPL_Pthread_LIBRARY_DIRS
    AND NOT TPL_Pthread_LIBRARIES
  )
  # Use CMake's Thread finder since it is a bit smarter in determining
  # whether pthreads is already built into the compiler and doesn't need
  # a library to link.
  find_package(Threads)
  #If Threads found a copy of pthreads make sure it is one of the cases the tribits
  #tpl system cannot handle.
  if (Threads_FOUND  AND  CMAKE_USE_PTHREADS_INIT)
    if(CMAKE_THREAD_LIBS_INIT STREQUAL "" OR  CMAKE_THREAD_LIBS_INIT  STREQUAL "-pthread")
      set(USE_THREADS TRUE)
      set(TPL_Pthread_INCLUDE_DIRS "")
      set(TPL_Pthread_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
      set(TPL_Pthread_LIBRARY_DIRS "")
    endif()
  endif()
endif()

if (USE_THREADS  AND  CMAKE_THREAD_LIBS_INIT  STREQUAL "")
  # Produce dummy Pthread::all_libs target and PthreadConfig.cmake file
  tribits_tpl_find_include_dirs_and_libraries(Pthread)
else()
  tribits_tpl_find_include_dirs_and_libraries( Pthread
    REQUIRED_HEADERS pthread.h
    REQUIRED_LIBS_NAMES pthread
    )
endif()
