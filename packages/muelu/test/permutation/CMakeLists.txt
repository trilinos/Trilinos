
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

# GH: remove all the permutation tests based on this?

TRIBITS_COPY_FILES_TO_BINARY_DIR(Permutation_cp
SOURCE_FILES test1.txt test2.txt test3.txt test4.txt test5.txt exp1.txt exp2.txt exp3.txt exp4.txt exp5.txt
)
