include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)

########## Probe for various thread configurations ##############

IF (CMAKE_HAVE_PTHREAD_H)

    # Probe for pthreads header file
    CHECK_INCLUDE_FILES("pthread.h" HAVE_TRIOS_PTHREAD_H)


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
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;return 0;}"
        HAVE_TRIOS_PTHREAD_MUTEX_T
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_mutex_init(&mutex,NULL);return 0;}"
        HAVE_TRIOS_PTHREAD_MUTEX_INIT
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_mutex_lock(&mutex);return 0;}"
        HAVE_TRIOS_PTHREAD_MUTEX_LOCK
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_mutex_unlock(&mutex);return 0;}"
        HAVE_TRIOS_PTHREAD_MUTEX_UNLOCK
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_mutex_destroy(&mutex);return 0;}"
        HAVE_TRIOS_PTHREAD_MUTEX_DESTROY
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_cond_t cond;return 0;}"
        HAVE_TRIOS_PTHREAD_COND_T
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_cond_t cond;pthread_cond_init(&cond,NULL);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_INIT
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_cond_t cond;pthread_cond_wait(&cond,&mutex);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_WAIT
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex;pthread_cond_t cond;struct timespec *abstime;pthread_cond_timedwait(&cond,&mutex,abstime);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_TIMEDWAIT
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_cond_t cond;pthread_cond_signal(&cond);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_SIGNAL
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_cond_t cond;pthread_cond_broadcast(&cond);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_BROADCAST
    )
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_cond_t cond;pthread_cond_destroy(&cond);return 0;}"
        HAVE_TRIOS_PTHREAD_COND_DESTROY
    )
    add_definitions(-D_GNU_SOURCE)  # needed to find RECURSIVE mutex initializer
    # Test for PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
    set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
    check_c_source_compiles(
        "#include <pthread.h>\nint main(){pthread_mutex_t mutex=PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;return 0;}"
        HAVE_TRIOS_PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP
    )

ENDIF()
