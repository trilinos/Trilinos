tribits_get_package_enable_status(Kokkos  KokkosEnable "")


macro(enable_warnings warnings)
    message(STATUS "Trilinos warnings enabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-W${warning} -Wno-error=${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


macro(disable_warnings warnings)
    message(STATUS "Trilinos warnings disabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-Wno-${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


macro(enable_errors errors)
    message(STATUS "Trilinos warnings-as-errors enabled: ${errors}")
    foreach(error ${errors})
        set(CMAKE_CXX_FLAGS "-Werror=${error} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()

message(STATUS "Adding '-std=c99' to C compiler flags for Zoltan")
set(Zoltan_C_FLAGS "-std=c99 ${Zoltan_C_FLAGS} ${CMAKE_C_FLAGS}")

IF (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  IF(WIN32)
    MESSAGE("-- " "Adding '/fp:precise' to C++ compiler flags because Trilinos needs it when using the Intel OneAPI C++ compiler.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:precise")
  ELSE()
    MESSAGE("-- " "Adding '-fp-model=precise' to C++ compiler flags because Trilinos needs it when using the Intel OneAPI C++ compiler.")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model=precise")
  ENDIF()
ENDIF()

IF (KokkosEnable)
  MESSAGE("-- " "Skip adding flags for OpenMP because Kokkos flags does that ...")
  SET(OpenMP_CXX_FLAGS_OVERRIDE " ")

  # There is a top-level CMAKE_CXX_FLAGS. Warnings and other flags get added
  # in the sub-scope of each individual package
  # Grab those variables here so we know what was the original top-level flags
  # and what are the CMAKE_CXX_FLAGS added afterwards for an individual package
  SET(TRILINOS_TOPLEVEL_CXX_FLAGS ${CMAKE_CXX_FLAGS})
  # NOTE: Setting TRILINOS_TOPLEVEL_CXX_FLAGS only has any impact if Kokkos is
  # being treated as an internal package.
ENDIF()

set(explicitly_disabled_warnings
    deprecated-declarations
    inline
)
set(upcoming_warnings
    aggressive-loop-optimizations
    array-bounds=2  # -Wall
    class-memaccess  # -Wall
    dangling-pointer=2  # -Wall
    # deprecated-copy  # -Wextra, lots of warnings
    implicit-fallthrough=3  # -Wextra
    maybe-uninitialized  # -Wall
    mismatched-new-delete  # -Wall
    pessimizing-move  # -Wall
    redundant-move  # -Wextra
    restrict
    #unused-parameter  # -Wextra, lots of warnings
    ${Trilinos_ADDITIONAL_WARNINGS}
)
set(promoted_warnings
    address
    aligned-new
    alloc-size  # -Wextra
    array-compare
    bool-compare
    bool-operation
    builtin-declaration-mismatch
    c++11-compat
    c++14-compat
    c++17compat
    c++20compat
    calloc-transposed-args
    cast-align
    cast-function-type  # -Wextra
    catch-value
    char-subscripts
    clobbered
    comment
    dangling-else
    dangling-reference  # -Wextra
    delete-non-virtual-dtor
    div-by-zero
    duplicate-decl-specifier
    empty-body
    enum-compare
    enum-int-mismatch
    expansion-to-defined  # -Wextra
    format
    format=1
    format-contains-nul
    format-diag
    format-extra-args
    format-overflow=1
    format-truncation=1
    format-zero-length
    frame-address
    ignored-qualifiers  # -Wextra
    implicit
    implicit-function-declaration
    implicit-int
    infinite-recursion
    init-self
    int-in-bool-context
    int-to-pointer-cast
    logical-not-parentheses
    main
    memset-elt-size
    memset-transposed-args
    misleading-indentation
    mismatched-dealloc
    missing-attributes
    missing-field-initializers  # -Wextra
    multistatement-macros
    narrowing
    nonnull
    nonnull-compare
    openmp-simd
    overloaded-virtual=1
    packed-not-aligned
    parentheses
    pointer-sign
    range-loop-construct
    reorder
    return-type
    self-move
    sequence-point
    shadow
    sign-compare
    sized-deallocation  # -Wextra
    sizeof-array-div
    sizeof-pointer-div
    sizeof-pointer-memaccess
    strict-aliasing
    strict-overflow=1
    string-compare  # -Wextra
    switch
    tautological-compare
    trigraphs
    type-limits
    uninitialized
    unknown-pragmas
    unused
    unused-but-set-parameter
    unused-but-set-variable
    unused-const-variable=1
    unused-function
    unused-label
    unused-local-typedefs
    unused-value
    unused-variable
    use-after-free=2
    variadic-macros
    vla-parameter
    volatile-register-var
    write-strings
    zero-length-bounds
)

include(CheckCXXCompilerFlag)

function(filter_valid_warnings warnings output)
    set(valid_warnings "")
    foreach(warning ${warnings})
        # Check if the compiler supports the warning flag
        string(CONCAT flag "-W" ${warning})
        check_cxx_compiler_flag("${flag}" COMPILER_SUPPORTS_${warning}_WARNING)

        if(COMPILER_SUPPORTS_${warning}_WARNING)
            list(APPEND valid_warnings "${warning}")
        endif()
    endforeach()
    set(${output} ${valid_warnings} PARENT_SCOPE)
endfunction()


function(filter_valid_warnings_as_errors warnings output)
    set(valid_warnings "")
    foreach(warning ${warnings})
        # Check if the compiler supports the warning-as-error flag
        string(CONCAT flag "-Werror=" ${warning})
        check_cxx_compiler_flag("${flag}" COMPILER_SUPPORTS_${warning}_WARNING_AS_ERROR)

        if(COMPILER_SUPPORTS_${warning}_WARNING_AS_ERROR)
            list(APPEND valid_warnings "${warning}")
        endif()
    endforeach()
    set(${output} ${valid_warnings} PARENT_SCOPE)
endfunction()

if("${Trilinos_WARNINGS_MODE}" STREQUAL "WARN")
    filter_valid_warnings("${upcoming_warnings}" upcoming_warnings)
    enable_warnings("${upcoming_warnings}")
    filter_valid_warnings_as_errors("${promoted_warnings}" promoted_warnings)
    enable_errors("${promoted_warnings}")
elseif("${Trilinos_WARNINGS_MODE}" STREQUAL "ERROR")
    filter_valid_warnings_as_errors("${promoted_warnings}" promoted_warnings)
    filter_valid_warnings_as_errors("${upcoming_warnings}" upcoming_warnings)
    enable_errors("${promoted_warnings};${upcoming_warnings}")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    disable_warnings("${explicitly_disabled_warnings}")
endif()


# Make sure that all deprecated packages are forcefully disabled (and that the variables controlling enablement are defined)
macro(force_disable_package packageName)
    tribits_filter_package_list_from_var(Trilinos_DEFINED_PACKAGES INTERNAL ON NONEMPTY packageSublist)
    set(Trilinos_ENABLE_${packageName} OFF CACHE BOOL "Enable ${packageName} (special setting for force-disable, should ALWAYS be `OFF`)")
    foreach(package ${packageSublist})
        set(${package}_ENABLE_${packageName} OFF CACHE BOOL "Enable ${packageName} support in ${package} (special setting for force-disable, should ALWAYS be `OFF`)")
    endforeach()
endmacro()

set(DEPRECATED_PACKAGES Amesos AztecOO Epetra EpetraExt Ifpack Intrepid Isorropia ML NewPackage Pliris PyTrilinos ShyLU_DDCore ThyraEpetraAdapters ThyraEpetraExtAdapters Triutils)
FOREACH(package ${DEPRECATED_PACKAGES})
  force_disable_package(${package})
ENDFOREACH()
