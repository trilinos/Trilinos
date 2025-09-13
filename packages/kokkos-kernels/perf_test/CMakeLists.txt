if (KokkosKernels_ENABLE_PERFTESTS)
    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

    KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../test_common)

    #At some point, we may wish to make this into real "performance
    #tests, " in the sense that they can pass or fail.  At that point, use
    #"CATEGORIES PERFORMANCE" to mark them as such.For now, we just
    #build the executables for manual use, but don't run the tests.  They
    #build correctly with or without MPI, but only run them with a single
    #MPI process.

    KOKKOSKERNELS_INCLUDE_DIRECTORIES(sparse)

    if(KokkosKernels_ENABLE_TESTS_AND_PERFSUITE)
        #Add RPS implementations of KK perf tests here
        KOKKOSKERNELS_ADD_EXECUTABLE(
            tracked_testing
            SOURCES KokkosKernelsTrackedTesting.cpp
            sparse/KokkosSparse_spmv_test.cpp
            blas/blas2/KokkosBlas2_gemv_tracked_perf_test.cpp
            blas/blas1/KokkosBlas_dot_tracked_perf_test.cpp
            blas/blas1/KokkosBlas_team_dot_tracked_perf_test.cpp
            blas/blas3/KokkosBlas3_gemm_tracked_perf_test.cpp
            PerfTestUtilities.cpp
            sparse/spmv/OpenMPSmartStatic_SPMV.cpp
            #sparse / KokkosSparse_spgemm_test.cpp
            )
    endif()

    ADD_COMPONENT_SUBDIRECTORY(batched)
    ADD_COMPONENT_SUBDIRECTORY(graph)
    ADD_COMPONENT_SUBDIRECTORY(sparse)
    ADD_COMPONENT_SUBDIRECTORY(blas)
    ADD_SUBDIRECTORY(performance)
    #ADD_SUBDIRECTORY(common)

endif()
