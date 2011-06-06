
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE TEUCHOS)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/../../../Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Teuchos_AbstractFactory.hpp
  Teuchos_AbstractFactoryStd.hpp
  Teuchos_Array.hpp
  Teuchos_ArrayRCP.hpp
  Teuchos_ArrayRCPDecl.hpp
  Teuchos_ArrayView.hpp
  Teuchos_ArrayViewDecl.hpp
  Teuchos_ArrayConversions.hpp
  Teuchos_Assert.hpp
  Teuchos_BLAS.hpp
  Teuchos_BLAS_types.hpp
  Teuchos_BLAS_wrappers.hpp
  Teuchos_CTimeMonitor.h
  Teuchos_Comm.hpp
  Teuchos_CommHelpers.hpp
  Teuchos_CommUtilities.hpp
  Teuchos_CommandLineProcessor.hpp
  Teuchos_CompObject.hpp
  Teuchos_CompileTimeAssert.hpp
  Teuchos_Condition.hpp
  Teuchos_ConfigDefs.hpp
  Teuchos_ConstNonconstObjectContainer.hpp
  Teuchos_ConstTypeTraits.hpp
  Teuchos_DataAccess.hpp
  Teuchos_DebugDefaultAsserts.hpp
  Teuchos_DefaultComm.hpp
  Teuchos_DefaultSerialComm.hpp
  Teuchos_Dependency.hpp
  Teuchos_DependencySheet.hpp
  Teuchos_Describable.hpp
  Teuchos_DummyObjectGetter.hpp
  Teuchos_ENull.hpp
  Teuchos_ErrorPolling.hpp
  Teuchos_Exceptions.hpp
  Teuchos_ExpandScalarTypeMacros.hpp
  Teuchos_ExplicitInstantiationHelpers.hpp
  Teuchos_F77_wrappers.h
  Teuchos_FancyOStream.hpp
  Teuchos_Flops.hpp
  Teuchos_FunctionObject.hpp
  Teuchos_GlobalMPISession.hpp
  Teuchos_Handle.hpp
  Teuchos_Handleable.hpp
  Teuchos_InvalidConditionException.hpp
  Teuchos_InvalidDependencyException.hpp
  Teuchos_LAPACK.hpp
  Teuchos_LAPACK_wrappers.hpp
  Teuchos_LabeledObject.hpp
  Teuchos_LocalTestingHelpers.hpp
  Teuchos_MPIComm.hpp
  Teuchos_MPIContainerComm.hpp
  Teuchos_MPISession.hpp
  Teuchos_MPITraits.hpp
  Teuchos_NullIteratorTraits.hpp
  Teuchos_Object.hpp
  Teuchos_ObjectBuilder.hpp
  Teuchos_OpaqueWrapper.hpp
  Teuchos_OrdinalTraits.hpp
  Teuchos_ParameterEntry.hpp
  Teuchos_ParameterEntryValidator.hpp
  Teuchos_ParameterList.hpp
  Teuchos_ParameterListAcceptor.hpp
  Teuchos_ParameterListAcceptorDefaultBase.hpp
  Teuchos_ParameterListAcceptorHelpers.hpp
  Teuchos_ParameterListExceptions.hpp
  Teuchos_ParameterListNonAcceptor.hpp
  Teuchos_PerformanceMonitorBase.hpp
  Teuchos_PerformanceMonitorUtils.hpp
  Teuchos_Polynomial.hpp
  Teuchos_PolynomialDecl.hpp
  Teuchos_PolynomialTraits.hpp
  Teuchos_PromotionTraits.hpp
  Teuchos_Ptr.hpp
  Teuchos_PtrDecl.hpp
  Teuchos_RCP.hpp
  Teuchos_RCPDecl.hpp
  Teuchos_RCPNode.hpp
  Teuchos_Range1D.hpp
  Teuchos_RawMPITraits.hpp
  Teuchos_ReductionOp.hpp
  Teuchos_ReductionOp.hpp
  Teuchos_ReductionOpHelpers.hpp
  Teuchos_RefCountPtr.hpp
  Teuchos_RefCountPtrDecl.hpp
  Teuchos_ScalarTraits.hpp
  Teuchos_ScalarTraitsCUDA.hpp
  Teuchos_ScalarTraitsDecl.hpp
  Teuchos_SerialDenseHelpers.hpp
  Teuchos_SerialDenseMatrix.hpp
  Teuchos_SerialDenseSolver.hpp
  Teuchos_SerialDenseVector.hpp
  Teuchos_SerialSpdDenseSolver.hpp
  Teuchos_SerialSymDenseMatrix.hpp
  Teuchos_SerializationTraits.hpp
  Teuchos_SerializationTraitsHelpers.hpp
  Teuchos_Serializer.hpp
  Teuchos_SerializerHelpers.hpp
  Teuchos_SimpleObjectTable.hpp
  Teuchos_StandardCatchMacros.hpp
  Teuchos_StandardCompositionMacros.hpp
  Teuchos_StandardConditions.hpp
  Teuchos_StandardDependencies.hpp
  Teuchos_StandardFunctionObjects.hpp
  Teuchos_StandardMemberCompositionMacros.hpp
  Teuchos_StandardParameterEntryValidators.hpp
  Teuchos_StaticSetupMacro.hpp
  Teuchos_StringToIntMap.hpp
  Teuchos_TableColumn.hpp
  Teuchos_TableEntry.hpp
  Teuchos_TableFormat.hpp
  Teuchos_TabularOutputter.hpp
  Teuchos_TestForException.hpp
  Teuchos_TestingHelpers.hpp
  Teuchos_Time.hpp
  Teuchos_TimeMonitor.hpp
  Teuchos_Tuple.hpp
  Teuchos_TypeNameTraits.hpp
  Teuchos_TypeTraits.hpp
  Teuchos_TwoDArray.hpp
  Teuchos_UnitTestBase.hpp
  Teuchos_UnitTestHarness.hpp
  Teuchos_UnitTestHelpers.hpp
  Teuchos_UnitTestRepository.hpp
  Teuchos_Utils.hpp
  Teuchos_ValidatorFactory.hpp
  Teuchos_VerboseObject.hpp
  Teuchos_VerboseObjectParameterListHelpers.hpp
  Teuchos_VerbosityLevel.hpp
  Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp
  Teuchos_Version.hpp
  Teuchos_Workspace.hpp
  Teuchos_any.hpp
  Teuchos_arrayArg.hpp
  Teuchos_as.hpp
  Teuchos_basic_oblackholestream.hpp
  Teuchos_dyn_cast.hpp
  Teuchos_getBaseObjVoidPtr.hpp
  Teuchos_getConst.hpp
  Teuchos_getRawPtr.hpp
  Teuchos_implicit_cast.hpp
  Teuchos_map.hpp
  Teuchos_oblackholestream.hpp
  Teuchos_stacktrace.hpp
  Teuchos_toString.hpp
  )

# Set sources

APPEND_SET(SOURCES
  Teuchos_BLAS.cpp
  Teuchos_CTimeMonitor.cpp
  Teuchos_CommandLineProcessor.cpp
  Teuchos_CompObject.cpp
  Teuchos_Describable.cpp
  Teuchos_Dependency.cpp
  Teuchos_DependencySheet.cpp
  Teuchos_ErrorPolling.cpp
  Teuchos_Flops.cpp
  Teuchos_GlobalMPISession.cpp
  Teuchos_LabeledObject.cpp
  Teuchos_MPIComm.cpp
  Teuchos_MPISession.cpp
  Teuchos_Object.cpp
  Teuchos_ParameterEntry.cpp
  Teuchos_ParameterList.cpp
  Teuchos_ParameterListAcceptorHelpers.cpp
  Teuchos_ParameterListAcceptor.cpp
  Teuchos_ParameterListAcceptorDefaultBase.cpp
  Teuchos_ParameterListNonAcceptor.cpp
  Teuchos_PerformanceMonitorUtils.cpp
  Teuchos_Ptr.cpp
  Teuchos_RCPNode.cpp
  Teuchos_Range1D.cpp
  Teuchos_ScalarTraits.cpp
  Teuchos_StandardConditions.cpp
  Teuchos_StandardDependencies.cpp
  Teuchos_StandardParameterEntryValidators.cpp
  Teuchos_StringToIntMap.cpp
  Teuchos_TableColumn.cpp
  Teuchos_TableEntry.cpp
  Teuchos_TableFormat.cpp
  Teuchos_TabularOutputter.cpp
  Teuchos_TestForException.cpp
  Teuchos_TestingHelpers.cpp
  Teuchos_Time.cpp
  Teuchos_TimeMonitor.cpp
  Teuchos_TypeNameTraits.cpp
  Teuchos_UnitTestBase.cpp
  Teuchos_UnitTestRepository.cpp
  Teuchos_Utils.cpp
  Teuchos_ValidatorFactory.cpp
  Teuchos_VerboseObject.cpp
  Teuchos_VerboseObjectParameterListHelpers.cpp
  Teuchos_VerbosityLevel.cpp
  Teuchos_VerbosityLevelCommandLineProcessorHelpers.cpp
  Teuchos_Workspace.cpp
  Teuchos_dyn_cast.cpp
  Teuchos_stacktrace.cpp
  )

# Set addition headers and source based on package-specific
# configuration options

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)
APPEND_SET(HEADERS
    Teuchos_ConditionXMLConverter.hpp
    Teuchos_ConditionXMLConverterDB.hpp
    Teuchos_DependencyXMLConverter.hpp
    Teuchos_DependencyXMLConverterDB.hpp
    Teuchos_ExpatHandlerAdapter.hpp
    Teuchos_FileInputSource.hpp
    Teuchos_FileInputStream.hpp
    Teuchos_FunctionObjectXMLConverter.hpp
    Teuchos_FunctionObjectXMLConverterDB.hpp
    Teuchos_HashSet.hpp
    Teuchos_Hashtable.hpp
    Teuchos_HashUtils.hpp
    Teuchos_Libxml2HandlerAdapter.hpp
    Teuchos_ParameterXMLFileReader.hpp
    Teuchos_ParameterEntryXMLConverter.hpp
    Teuchos_ParameterEntryXMLConverterDB.hpp
    Teuchos_StandardConditionXMLConverters.hpp
    Teuchos_StandardDependencyXMLConverters.hpp
    Teuchos_StandardFunctionObjectXMLConverters.hpp
    Teuchos_StandardParameterEntryXMLConverters.hpp
    Teuchos_StandardValidatorXMLConverters.hpp
    Teuchos_StringInputSource.hpp
    Teuchos_StringInputStream.hpp
    Teuchos_StrUtils.hpp
    Teuchos_TreeBuildingXMLHandler.hpp
    Teuchos_ValidatorMaps.hpp
    Teuchos_ValidatorXMLConverter.hpp
    Teuchos_ValidatorXMLConverterDB.hpp
    Teuchos_XMLConditionExceptions.hpp
    Teuchos_XMLDependencyExceptions.hpp
    Teuchos_XMLFunctionObjectExceptions.hpp
    Teuchos_XMLInputSource.hpp
    Teuchos_XMLInputStream.hpp
    Teuchos_XMLObject.hpp
    Teuchos_XMLObjectImplem.hpp
    Teuchos_XMLParameterListExceptions.hpp
    Teuchos_XMLParameterListHelpers.hpp
    Teuchos_XMLParameterListReader.hpp
    Teuchos_XMLParameterListWriter.hpp
    Teuchos_XMLParser.hpp
    )

  APPEND_SET(SOURCES
    Teuchos_ConditionXMLConverter.cpp
    Teuchos_ConditionXMLConverterDB.cpp
    Teuchos_DependencyXMLConverter.cpp
    Teuchos_DependencyXMLConverterDB.cpp
    Teuchos_ExpatHandlerAdapter.cpp
    Teuchos_FileInputSource.cpp
    Teuchos_FileInputStream.cpp
    Teuchos_FunctionObjectXMLConverter.cpp
    Teuchos_FunctionObjectXMLConverterDB.cpp
    Teuchos_HashUtils.cpp
    Teuchos_Libxml2HandlerAdapter.cpp
    Teuchos_ParameterXMLFileReader.cpp
    Teuchos_ParameterEntryXMLConverter.cpp
    Teuchos_ParameterEntryXMLConverterDB.cpp
    Teuchos_StandardDependencyXMLConverters.cpp
    Teuchos_StandardParameterEntryXMLConverters.cpp
    Teuchos_StandardValidatorXMLConverters.cpp
    Teuchos_StandardConditionXMLConverters.cpp
    Teuchos_StrUtils.cpp
    Teuchos_StringInputSource.cpp
    Teuchos_StringInputStream.cpp
    Teuchos_TreeBuildingXMLHandler.cpp
    Teuchos_ValidatorMaps.cpp
    Teuchos_ValidatorXMLConverter.cpp
    Teuchos_ValidatorXMLConverterDB.cpp
    Teuchos_XMLInputSource.cpp
    Teuchos_XMLInputStream.cpp
    Teuchos_XMLObject.cpp
    Teuchos_XMLObjectImplem.cpp
    Teuchos_XMLParameterListHelpers.cpp
    Teuchos_XMLParameterListReader.cpp
    Teuchos_XMLParameterListWriter.cpp
    Teuchos_XMLParser.cpp
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_MPI)

  APPEND_SET(HEADERS
    Teuchos_DefaultMpiComm.hpp
    Teuchos_MpiReductionOpSetter.hpp
    )

  APPEND_SET(SOURCES
    Teuchos_MpiReductionOpSetter.cpp
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Boost)

  APPEND_SET(HEADERS
    Teuchos_RCPBoostSharedPtrConversions.hpp
    Teuchos_RCPBoostSharedPtrConversionsDecl.hpp
    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  # Must suppress strick aliasing warnings for this file!
  SET_SOURCE_FILES_PROPERTIES(
    Teuchos_stacktrace.cpp
    PROPERTIES COMPILE_FLAGS -fno-strict-aliasing
    )
ENDIF()

PACKAGE_ADD_LIBRARY(
  teuchos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOS_LIB_EXPORTS_MODE
  )
