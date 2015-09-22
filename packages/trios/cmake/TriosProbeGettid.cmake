
########## CHECK FOR HEADER FILES ############

INCLUDE(CheckIncludeFiles)

# Probe for syscall header files
CHECK_INCLUDE_FILES("unistd.h" HAVE_TRIOS_UNISTD_H)
CHECK_INCLUDE_FILES("syscall.h" HAVE_TRIOS_SYSCALL_H)

########## CHECK FOR FUNCTIONS ############

INCLUDE(CheckSymbolExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckCSourceCompiles)

# syscall()
IF (HAVE_TRIOS_UNISTD_H)
    CHECK_FUNCTION_EXISTS(syscall HAVE_TRIOS_SYSCALL)
ENDIF (HAVE_TRIOS_UNISTD_H)

# SYS_gettid
IF (HAVE_TRIOS_SYSCALL_H)
    CHECK_SYMBOL_EXISTS(SYS_gettid "syscall.h" HAVE_TRIOS_SYS_GETTID)
ENDIF (HAVE_TRIOS_SYSCALL_H)

IF (HAVE_TRIOS_SYSCALL)
    IF (HAVE_TRIOS_SYS_GETTID)
        check_c_source_compiles(
            "#include <unistd.h>\n#include <syscall.h>\nint main(){syscall(SYS_gettid);return 0;}"
            HAVE_TRIOS_GETTID
        )
    ENDIF (HAVE_TRIOS_SYS_GETTID)
ENDIF (HAVE_TRIOS_SYSCALL)
