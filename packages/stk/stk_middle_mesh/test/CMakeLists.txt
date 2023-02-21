enable_testing()

find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB HELPER_LIB_SRCS util/*.cc)
add_library(test_helper_lib ${HELPER_LIB_SRCS})
target_include_directories(test_helper_lib PRIVATE ${GTEST_INCLUDE_DIRS} "${CMAKE_CURRENT_SOURCE_DIR}/utils")
target_link_libraries(test_helper_lib ${GTEST_BOTH_LIBRARIES} triangle_lib)

file(GLOB SRCS *.cc)

add_executable(unit_tests ${SRCS})
target_include_directories(unit_tests PRIVATE ${GTEST_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(unit_tests ${GTEST_BOTH_LIBRARIES} triangle_lib test_helper_lib)



file(GLOB VERIFICATION_SRCS verification_tests/*.cc)

message("verification test sources = ${VERIFICATION_SRCS}")
add_executable(verification_tests ${VERIFICATION_SRCS} gtest_main.cc) 
target_include_directories(verification_tests PRIVATE ${GTEST_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(verification_tests ${GTEST_BOTH_LIBRARIES} triangle_lib test_helper_lib)


#add_test(AllUnitsTests unit_tests)
add_test(NAME UnitTests_np1 COMMAND unit_tests)
add_test(NAME UnitTests_np2 COMMAND mpiexec -n 2 unit_tests)
add_test(NAME UnitTests_np3 COMMAND mpiexec -n 3 unit_tests)
add_test(NAME UnitTests_np4 COMMAND mpiexec -n 4 unit_tests)
add_test(NAME UnitTests_np5 COMMAND mpiexec -n 5 unit_tests)
add_test(NAME UnitTests_np6 COMMAND mpiexec -n 6 unit_tests)
add_test(NAME UnitTests_np7 COMMAND mpiexec -n 7 unit_tests)
add_test(NAME UnitTests_np8 COMMAND mpiexec -n 8 unit_tests)
