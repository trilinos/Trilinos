find_package(sierra-cmake-utils REQUIRED)


function(add_custom_library product lib_name src_glob header_glob include_dir)
    if (lib_name MATCHES "^INTERFACE_")
        set(lib_type "INTERFACE")
        string(REPLACE "INTERFACE_" "" lib_name ${lib_name})
        add_library(${lib_name} INTERFACE)
    else()
        set(lib_type "PUBLIC")
        add_library(${lib_name})
        file(GLOB ${lib_name}_sources CONFIGURE_DEPENDS ${src_glob})
        target_sources(${lib_name} PRIVATE ${${lib_name}_sources})
    endif()

    if (ARGN)
        foreach(lib IN LISTS ARGN)
            if (lib MATCHES "^PRIVATE_.*")
                string(REPLACE "PRIVATE_" "" private_lib ${lib})
                target_link_libraries(${lib_name} PRIVATE ${private_lib})
            elseif (lib MATCHES "^PUBLIC_.*")
                string(REPLACE "PUBLIC_" "" public_lib ${lib})
                target_link_libraries(${lib_name} PUBLIC ${public_lib})
            elseif (lib MATCHES "^INTERFACE_.*")
                string(REPLACE "INTERFACE_" "" interface_lib ${lib})
                target_link_libraries(${lib_name} INTERFACE ${interface_lib})
            endif()
        endforeach()
    endif()

    file(GLOB ${lib_name}_headers CONFIGURE_DEPENDS ${header_glob})

    # Only advertise the per-library include directory when the library
    # actually ships headers.  A source-free umbrella library (e.g. krino_lib)
    # installs no headers, so its include/${product}/${include_dir} path is
    # never created; exporting it would make downstream find_package() fail
    # with "Imported target includes non-existent path".  Such umbrella
    # libraries re-export the include dirs of the real libraries they link
    # transitively, so they need no include directory of their own.
    if (${lib_name}_headers)
        target_include_directories(${lib_name} ${lib_type}
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${include_dir}>
            $<INSTALL_INTERFACE:include/${product}/${include_dir}>
        )
    endif()

    target_sources(${lib_name} ${lib_type}
        FILE_SET ${lib_name}_headers
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
        FILES ${${lib_name}_headers}
    )

    install(
        TARGETS ${lib_name}
        EXPORT ${product}Targets
        FILE_SET ${lib_name}_headers
        DESTINATION include/${product}
        INCLUDES DESTINATION include/${product}
    )
endfunction()


function(add_custom_library_in_dir product lib_name dir_name)
    add_custom_library(
        ${product}
        ${lib_name}
        "${dir_name}/*.cpp"
        "${dir_name}/*.hpp"
        "${dir_name}"
         ${ARGN}
    )
endfunction()


function(add_custom_executable exe_name src_glob)
    add_executable(${exe_name})

    file(GLOB ${exe_name}_sources CONFIGURE_DEPENDS ${src_glob})
    target_sources(${exe_name} PRIVATE ${${exe_name}_sources})

    if (ARGN)
        foreach(arg IN LISTS ARGN)
            if (arg MATCHES "^PUBLIC_.*")
                string(REPLACE "PUBLIC_" "" public_lib ${arg})
                target_link_libraries(${exe_name} PUBLIC ${public_lib})
            elseif (arg MATCHES "^PRIVATE_.*")
                string(REPLACE "PRIVATE_" "" private_lib ${arg})
                target_link_libraries(${exe_name} PRIVATE ${private_lib})
            elseif (arg MATCHES "^INCLUDE_.*")
                string(REPLACE "INCLUDE_" "" include_dir ${arg})
                target_include_directories(${exe_name} PRIVATE ${include_dir})
            endif()
        endforeach()
    endif()

    install(TARGETS ${exe_name})
endfunction()


function(add_unit_test_on_num_procs exe_name num_procs)
    sierra_add_test(NAME ${exe_name} EXECUTABLE ${exe_name} RANKS ${num_procs})
endfunction()


function(add_parallel_unit_test_executable exe_name src_glob)
    if (BUILD_TESTS)
        add_custom_executable(${exe_name} ${src_glob} ${ARGN})
        add_unit_test_on_num_procs(${exe_name} 1)
        add_unit_test_on_num_procs(${exe_name} 2)
        add_unit_test_on_num_procs(${exe_name} 3)
        add_unit_test_on_num_procs(${exe_name} 4)
    endif()
endfunction()


function(add_unit_test_executable exe_name src_glob)
    if (BUILD_TESTS)
        add_custom_executable(${exe_name} ${src_glob} ${ARGN})
        sierra_add_test(NAME ${exe_name} EXECUTABLE ${exe_name})
    endif()
endfunction()


function(add_regression_test_executable exe_name src_glob)
    if (BUILD_TESTS)
        add_custom_executable(${exe_name} ${src_glob} ${ARGN})
    endif()
endfunction()


function(add_custom_executable_with_build_info prefix exe_name src_glob)
    add_custom_executable(${exe_name} ${src_glob} ${ARGN})
    
    string(TIMESTAMP build_date)
    cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
     
    target_compile_definitions(${exe_name} PRIVATE
        ${prefix}_BUILD_DATE="${build_date}"
        ${prefix}_BUILD_MACH="${hostname}"
        ${prefix}_BUILD_OS1="${CMAKE_HOST_SYSTEM_NAME}"
        ${prefix}_BUILD_OS2="${CMAKE_HOST_SYSTEM_VERSION}"
        ${prefix}_BUILD_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}-$<CONFIG>"
    )
endfunction()

