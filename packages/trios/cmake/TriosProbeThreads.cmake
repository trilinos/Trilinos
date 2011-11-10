include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)

########## Probe for various thread configurations ##############

add_definitions(-D_GNU_SOURCE)  # needed to find RECURSIVE mutex initializer
# Test for PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
check_c_source_compiles(
    "#include <pthread.h>\nint main(){pthread_mutex_t mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;return 0;}"
    HAVE_TRIOS_PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
)
# Test for a way to yield besides usleep(0)
check_c_source_compiles(
    "#include <pthread.h>\nint main(){pthread_yield();return 0;}"
    HAVE_TRIOS_PTHREAD_YIELD
)
check_c_source_compiles(
    "#include <pthread.h>\nint main(){pthread_yield_np();return 0;}"
    HAVE_TRIOS_PTHREAD_YIELD_NP
)
check_c_source_compiles(
    "#include <sched.h>\nint main(){sched_yield();return 0;}"
    HAVE_TRIOS_SCHED_YIELD
)
