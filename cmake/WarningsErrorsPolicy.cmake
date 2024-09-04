set(warnings shadow parentheses sign-compare unused-variable)


macro(disable_warnings_for_deprecated_packages)
    message(STATUS "Disabling all warnings/errors for deprecated packages")
    foreach(package deprecated_packages)
        set(${package}_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endforeach()
endmacro()


if("${Trilinos_WARNINGS}" STREQUAL "WARN")
    message(STATUS "Trilinos warnings enabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-W${warning} -Wno-error=${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
    disable_warnings_for_deprecated_packages()
elseif("${Trilinos_WARNINGS}" STREQUAL "ERROR")
    message(STATUS "Trilinos warnings-as-errors enabled: ${warnings}")
    foreach(warning ${warnings})
        set(CMAKE_CXX_FLAGS "-W${warning} -Wno-error=${warning} ${CMAKE_CXX_FLAGS}")
    endforeach()
    disable_warnings_for_deprecated_packages()
endif()
