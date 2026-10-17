IF(COMMAND TRIBITS_PACKAGE_DECL)
    message("Building PyKrino as part of Krino Trilinos package")
    # -----------------------------------------------------------------------
    # In-tree build under the TriBITS/Trilinos-package path (add_subdirectory
    # from krino/CMakeLists.txt). Original behavior, unchanged.
    # -----------------------------------------------------------------------
    # Get the python version
    EXECUTE_PROCESS(COMMAND ${Python3_EXECUTABLE} -c
                            "import sys; print(sys.version_info.major)"
      OUTPUT_VARIABLE PYTHON_MAJOR_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    EXECUTE_PROCESS(COMMAND ${Python3_EXECUTABLE} -c
                            "import sys; print(sys.version_info.minor)"
      OUTPUT_VARIABLE PYTHON_MINOR_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    SET(PYTHON_VERSION ${PYTHON_MAJOR_VERSION}.${PYTHON_MINOR_VERSION})

    SET(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})

    SET(PyKrino_INSTALL_DIR
        lib/python${PYTHON_VERSION}/site-packages/pykrino
        CACHE STRING "The path where PyKrino will be installed"
        )

    SET(PyKrino_PYTHONPATH "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/src")

    ADD_SUBDIRECTORY(src)

    TRIBITS_ADD_TEST_DIRECTORIES(test)
    FILE(COPY test DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else()
    MESSAGE("*** Building py-krino as a stand-alone Spack package. ***")
    include(SierraCmakeCode.cmake)
endif()
