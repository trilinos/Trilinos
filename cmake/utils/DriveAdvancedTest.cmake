INCLUDE(PrintVar)

FUNCTION(DRIVE_ADVANCED_TEST)

  SET(TEST_SEP
    "================================================================================")

  SET(OUTPUT_SEP
    "--------------------------------------------------------------------------------")

  SET(TEST_PASSED TRUE)
  
  MATH(EXPR LAST_CMND_IDX ${NUM_CMNDS}-1)
  
  FOREACH ( CMND_IDX RANGE ${LAST_CMND_IDX})

    MESSAGE("\n${TEST_SEP}\n")
    MESSAGE("TEST_${CMND_IDX}\n")

    MESSAGE("Runing: ${TEST_CMND_${CMND_IDX}} ...\n")
    SET(EXEC_CMND COMMAND ${TEST_CMND_${CMND_IDX}})

    IF (TEST_OUTPUT_FILE_${CMND_IDX})
      MESSAGE("  Writing output to file ${TEST_OUTPUT_FILE_${CMND_IDX}} ...")
      SET(EXEC_CMND ${EXEC_CMND} OUTPUT_FILE ${TEST_OUTPUT_FILE_${CMND_IDX}})
    ENDIF()

    EXECUTE_PROCESS(${EXEC_CMND} RESULT_VARIABLE EXEC_RESULT)

    MESSAGE("\n${OUTPUT_SEP}\n")

    IF (EXEC_RESULT EQUAL 0)
      MESSAGE("TEST_${CMND_IDX} PASSED: Return code = 0\n")
    ELSE()
      MESSAGE("TEST_${CMND_IDX} FAILED: Return code = ${EXEC_RESULT}\n")
      SET(TEST_PASSED FALSE)
    ENDIF()
  
  ENDFOREACH()
        
  MESSAGE("\n${TEST_SEP}\n")
    IF (TEST_PASSED)
      MESSAGE("FINAL RESULT: TEST PASSED")
    ELSE()
      MESSAGE("FINAL RESULT: TEST FAILED")
    ENDIF()
        
  MESSAGE("\n${TEST_SEP}\n")

ENDFUNCTION()
