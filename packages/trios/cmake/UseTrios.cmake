#
#
# See the file TriosConfig.cmake for all variables set when 
# executing the FindPackage(Trios) command.  
# 
# Compilers
#   Trios_CXX_COMPILER
#   Trios_C_COMPILER
#   Trios_FORTRAN_COMPILER
#  
# Compiler Flags
#   Trios_CXX_FLAGS
#   Trios_C_FLAGS
#   Trios_FORTRAN_FLAGS
#   Trios_EXTRA_LD_FLAGS
# 
# Paths
#   Trios_INCLUDE_DIRS
#   Trios_LIBRARY_DIRS
#   Trios_LIBRARIES
#   Trios_TPL_INCLUDE_DIRS
#   Trios_TPL_LIBRARY_DIRS
#   Trios_TPL_LIBRARIES
#
IF(NOT Trios_USE_FILE_INCLUDED)
    SET(Trios_USE_FILE_INCLUDED 1)
    
    if (Trios_FOUND)
    
        # Add compiler flags 
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${Trios_C_FLAGS}")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Trios_CXX_FLAGS}")
        SET(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} ${Trios_FORTRAN_FLAGS}")

        # Add include directories 
        INCLUDE_DIRECTORIES(${Trios_INCLUDE_DIRS})
        INCLUDE_DIRECTORIES(${Trios_TPL_INCLUDE_DIRS})

        # Add link directories
        LINK_DIRECTORIES(${Trios_LIBRARY_DIRS})
        LINK_DIRECTORIES(${Trios_TPL_LIBRARY_DIRS})

    
        # Add libraries directories
        SET(CMAKE_EXTRA_LIBS "${CMAKE_EXTRA_LIBS}")

        # Add path to CMAKE files
        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${Trios_DIR}")
    
    ENDIF (Trios_FOUND)

ENDIF(NOT Trios_USE_FILE_INCLUDED)
