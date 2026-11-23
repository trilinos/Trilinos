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

# Windows toolchains with the MSVC C++ ABI need /EHsc so throw/catch work in all
# packages and tests.  Covers Visual Studio and Ninja generators, cl.exe, clang-cl,
# Intel (classic icpc), and IntelLLVM on Windows - not MinGW/Cygwin (different models).
IF(WIN32 AND NOT CYGWIN AND ${PROJECT_NAME}_ENABLE_CXX)
  SET(_TRILINOS_WINDOWS_MSVC_ABI_CXX FALSE)
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    SET(_TRILINOS_WINDOWS_MSVC_ABI_CXX TRUE)
  ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND "x${CMAKE_CXX_SIMULATE_ID}" STREQUAL "xMSVC")
    SET(_TRILINOS_WINDOWS_MSVC_ABI_CXX TRUE)
  ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    SET(_TRILINOS_WINDOWS_MSVC_ABI_CXX TRUE)
  ELSEIF(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    SET(_TRILINOS_WINDOWS_MSVC_ABI_CXX TRUE)
  ENDIF()
  IF(_TRILINOS_WINDOWS_MSVC_ABI_CXX)
    IF(CMAKE_CXX_FLAGS MATCHES "(^| )/EHa( |$)")
      MESSAGE(WARNING
        "CMAKE_CXX_FLAGS contains /EHa; Trilinos will not add /EHsc (exception model left to user).")
    ELSEIF(CMAKE_CXX_FLAGS MATCHES "(^| )/EHsc-( |$)")
      MESSAGE(WARNING
        "CMAKE_CXX_FLAGS contains /EHsc-; Trilinos will not add /EHsc (exceptions explicitly disabled).")
    ELSEIF(NOT CMAKE_CXX_FLAGS MATCHES "(^| )/EHsc( |$)")
      MESSAGE("-- "
        "Adding '/EHsc' for C++ exception handling on Windows (${CMAKE_CXX_COMPILER_ID})")
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
      ADD_COMPILE_OPTIONS($<$<COMPILE_LANGUAGE:CXX>:/EHsc>)
    ELSE()
      # CMAKE_CXX_FLAGS already has /EHsc (e.g. CMake MSVC defaults); still set
      # directory compile options so a later -DCMAKE_CXX_FLAGS= override cannot
      # drop exceptions for targets in this project.
      ADD_COMPILE_OPTIONS($<$<COMPILE_LANGUAGE:CXX>:/EHsc>)
    ENDIF()
  ENDIF()
  UNSET(_TRILINOS_WINDOWS_MSVC_ABI_CXX)
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

# Listing of warning options for GCC
# https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html

set(explicitly_disabled_warnings
    deprecated-declarations
    inline
)
set(upcoming_warnings
    # Missing for -Wall
    array-bounds=1
    array-parameter=2
    dangling-pointer=2
    maybe-uninitialized
    pessimizing-move
    restrict

    # Missing for -Wextra
    # deprecated-copy  # lots of warnings
    implicit-fallthrough=3
    # unused-parameter  # lots of warnings
    ${Trilinos_ADDITIONAL_WARNINGS}
)
set(promoted_warnings

    # All the warnings set by -Wall
    address
    aligned-new
    # array-bounds=1
    array-compare
    # array-parameter=2
    bool-compare
    bool-operation
    c++11-compat
    c++14-compat
    c++17compat
    c++20compat
    catch-value
    char-subscripts
    class-memaccess
    comment
    dangling-else
    # dangling-pointer=2
    delete-non-virtual-dtor
    duplicate-decl-specifier
    enum-compare
    enum-int-mismatch
    format=1
    format-contains-nul
    format-diag
    format-extra-args
    format-overflow=1
    format-truncation=1
    format-zero-length
    frame-address
    implicit
    implicit-function-declaration
    implicit-int
    infinite-recursion
    init-self
    int-in-bool-context
    logical-not-parentheses
    main
    # maybe-uninitialized
    memset-elt-size
    memset-transposed-args
    misleading-indentation
    mismatched-dealloc
    mismatched-new-delete
    missing-attributes
    # missing-braces  # only turned on for C/ObjC by -Wall
    multistatement-macros
    narrowing
    nonnull
    nonnull-compare
    openmp-simd
    overloaded-virtual=1
    packed-not-aligned
    parentheses
    # pessimizing-move
    pointer-sign
    range-loop-construct
    reorder
    # restrict
    return-type
    self-move
    sequence-point
    sign-compare
    sizeof-array-div
    sizeof-pointer-div
    sizeof-pointer-memaccess
    strict-aliasing
    switch
    tautological-compare
    trigraphs
    uninitialized
    unknown-pragmas
    unused
    unused-but-set-variable
    unused-const-variable=1
    unused-function
    unused-label
    unused-local-typedefs
    unused-value
    unused-variable
    use-after-free=2
    vla-parameter
    volatile-register-var
    zero-length-bounds

    # All the warnings set by -Wextra
    # absolute-value
    alloc-size
    calloc-transposed-args
    cast-function-type
    clobbered
    dangling-reference
    # deprecated-copy
    empty-body
    # enum-conversion
    expansion-to-defined
    ignored-qualifiers
    # implicit-fallthrough=3
    # maybe-uninitialized is already in -Wall
    missing-field-initializers
    # missing-parameter-name
    # missing-parameter-type
    # old-style-declaration
    # multiple-parameter-fwd-decl-lists
    # override-init
    redundant-move
    # shift-negative-value (in C++11 to C++17 and in C99 and newer)
    # sign-compare is already in -Wall
    sized-deallocation
    string-compare
    type-limits
    # uninitialized is already in -Wall
    # unterminated-string-initialization
    # unused-parameter
    unused-but-set-parameter

    # other warning flags
    aggressive-loop-optimizations
    builtin-declaration-mismatch
    cast-align
    div-by-zero
    format
    int-to-pointer-cast
    shadow
    strict-overflow=1
    variadic-macros
    write-strings
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
