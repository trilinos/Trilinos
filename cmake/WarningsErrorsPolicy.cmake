macro(disable_warnings_for_deprecated_packages)
    message(STATUS "Disabling all warnings/errors for deprecated packages")
    foreach(package ${deprecated_packages})
        set(${package}_CXX_FLAGS "-w ${${package}_CXX_FLAGS}")
    endforeach()
endmacro()


macro(enable_warnings warnings)
    message(STATUS "Trilinos warnings enabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-W${warning} -Wno-error=${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


macro(enable_errors errors)
    message(STATUS "Trilinos warnings-as-errors enabled: ${errors}")
    foreach(error ${errors})
        set(CMAKE_CXX_FLAGS "-Werror=${error} ${CMAKE_CXX_FLAGS}")
    endforeach()
endmacro()


set(upcoming_warnings shadow ${Trilinos_ADDITIONAL_WARNINGS})
set(promoted_warnings parentheses sign-compare unused-variable)


if("${Trilinos_WARNINGS_MODE}" STREQUAL "WARN")
    enable_warnings("${upcoming_warnings}")
    enable_errors("${promoted_warnings}")
    disable_warnings_for_deprecated_packages()
elseif("${Trilinos_WARNINGS_MODE}" STREQUAL "ERROR")
    enable_errors("${promoted_warnings};${upcoming_warnings}")
    disable_warnings_for_deprecated_packages()
endif()
