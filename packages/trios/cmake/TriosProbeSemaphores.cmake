########## CHECK FOR HEADER FILES ############

INCLUDE(CheckIncludeFiles)

# Probe for semaphore header file
CHECK_INCLUDE_FILES("semaphore.h" HAVE_TRIOS_SEMAPHORE_H)

########## CHECK FOR FUNCTIONS ############

include(CheckCSourceCompiles)
INCLUDE(CheckCSourceRuns)


IF (HAVE_TRIOS_SEMAPHORE_H)

    check_c_source_compiles(
        "#include <semaphore.h>\nint main(){sem_t lock;return 0;}"
        HAVE_TRIOS_SEM_T
    )

    # Probe for unnamed semaphore implementation
        SET(SOURCE
        "
        #include <fcntl.h>
        #include <semaphore.h>
        int main()
        {
            sem_t lock;
            if (sem_init(&lock, 0, 1) == -1) {
                return(1);
            }
            if (sem_wait(&lock) == -1) {
                return(2);
            }
            if (sem_post(&lock) == -1) {
                return(3);
            }
            if (sem_destroy(&lock) == -1) {
                return(4);
            }

            return 0;
        }
        "
        )

        CHECK_C_SOURCE_RUNS("${SOURCE}" HAVE_TRIOS_UNNAMED_SEMAPHORES)

    # Probe for named semaphore implementation
        SET(SOURCE
        "
        #include <fcntl.h>
        #include <semaphore.h>
        int main()
        {
            sem_t *lock;
            lock=sem_open(\"/trios.sem_test\", O_CREAT, 0600, 1);
            if (lock == SEM_FAILED) {
                return(1);
            }
            if (sem_wait(lock) == -1) {
                return(2);
            }
            if (sem_post(lock) == -1) {
                return(3);
            }
            if (sem_close(lock) == -1) {
                return(4);
            }
            if (sem_unlink(\"/trios.sem_test\") == -1) {
                return(5);
            }

            return 0;
        }
        "
        )

        CHECK_C_SOURCE_RUNS("${SOURCE}" HAVE_TRIOS_NAMED_SEMAPHORES)

ENDIF(HAVE_TRIOS_SEMAPHORE_H)
