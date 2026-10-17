# ---------------------------------------------------------------------------
# Krino-specific wrappers around commonFunctions.cmake
# ---------------------------------------------------------------------------

function(krino_apply_common_compile_definitions target_name)
    target_compile_definitions(${target_name} PUBLIC KRINO_BUILT_IN_SIERRA)

    if ("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
        target_compile_definitions(${target_name} PUBLIC Build64)
    endif()
endfunction()

function(krino_apply_library_warnings target_name)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        target_compile_options(${target_name} PUBLIC
            $<$<COMPILE_LANGUAGE:C>:-Wshadow -Winconsistent-missing-override>
        )
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_compile_options(${target_name} PUBLIC
            $<$<COMPILE_LANGUAGE:C>:-Wshadow>
        )
    endif()
endfunction()

function(krino_apply_root_include_dirs target_name visibility)
    target_include_directories(${target_name} ${visibility}
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include/krino>
    )
endfunction()

function(krino_add_library lib_name dir_name)
    add_custom_library_in_dir(
        krino
        ${lib_name}
        ${dir_name}
        ${ARGN}
    )

    krino_apply_root_include_dirs(${lib_name} PUBLIC)
    krino_apply_common_compile_definitions(${lib_name})
    krino_apply_library_warnings(${lib_name})
endfunction()

function(krino_add_interface_library lib_name dir_name)
    add_custom_library(
        krino
        INTERFACE_${lib_name}
        ""
        "${dir_name}/*.hpp"
        "${dir_name}"
        ${ARGN}
    )

    krino_apply_root_include_dirs(${lib_name} INTERFACE)
endfunction()

function(krino_add_executable exe_name src_glob)
    add_custom_executable(
        ${exe_name}
        "${src_glob}"
        ${ARGN}
    )

    krino_apply_root_include_dirs(${exe_name} PUBLIC)
    krino_apply_common_compile_definitions(${exe_name})
endfunction()
