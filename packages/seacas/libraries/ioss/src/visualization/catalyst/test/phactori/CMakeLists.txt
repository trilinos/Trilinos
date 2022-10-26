
function(generatePhactoriUnitTests)
  set(PhactoriClassUnitTests
    test_PhactoriClipPlaneOperation
    test_PhactoriContourOperation
    test_PhactoriDataArtifactMetaDataControl
    test_PhactoriGeometricCellSampler1
    test_PhactoriMarkCellSurfaceStatus2
    test_PhactoriPointSourceFromJsonList
    test_PhactoriSampledCellInfo
    test_PhactoriSegmentCellSampler3
    test_PhactoriSegment
    test_PhactoriSliceWithPlaneOperation
    test_PhactoriThresholdOperation
    test_PhactoriVtkCellOperations
    test_PhactoriVtkDataExportOperation
  )

  set(PhactoriTestDirectory
    ${CMAKE_CURRENT_BINARY_DIR}/PhactoriTests)

  set(PhactoriUnitTestsDirectory
    ${PhactoriTestDirectory}/UnitTests)

  file(MAKE_DIRECTORY ${PhactoriUnitTestsDirectory})

  createPhactoriUnitTests("${PhactoriClassUnitTests}"
    ${PhactoriUnitTestsDirectory})
endfunction()

function(createPhactoriUnitTests Modules PhctrUnitTestDir)
  set (CURRENT_PYTHON_PATH $ENV{PYTHONPATH})
  ####set (LOCAL_PARAVIEW_VAR $ENV{PARAVIEW_ROOT})
  set (LOCAL_PARAVIEW_VAR ${PARAVIEW_INSTALL_DIR})
  set (PHACTORI_DIR ${CMAKE_CURRENT_BINARY_DIR}/../../phactori)
  set (PHACTORI_OPERATION_DIR ${PHACTORI_DIR}/Operation)
  set (PVBATCH_EXECUTABLE "${LOCAL_PARAVIEW_VAR}/bin/pvbatch")
  foreach(OnePhactoriClassTest ${Modules})
    add_test(NAME Phactori-Unit-${OnePhactoriClassTest}
      COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PHACTORI_DIR}:${CURRENT_PYTHON_PATH}
        ${PVBATCH_EXECUTABLE} ${PHACTORI_OPERATION_DIR}/${OnePhactoriClassTest}.py
        WORKING_DIRECTORY ${PhctrUnitTestDir}
    )
  endforeach()
endfunction()

generatePhactoriUnitTests()
