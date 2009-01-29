
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${PROJECT_SOURCE_DIR}/Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
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
  Teuchos_ConfigDefs.hpp
  Teuchos_ConstNonconstObjectContainer.hpp
  Teuchos_ConstTypeTraits.hpp
  Teuchos_DataAccess.hpp
  Teuchos_DefaultComm.hpp
  Teuchos_DefaultSerialComm.hpp
  Teuchos_Describable.hpp
  Teuchos_ENull.hpp
  Teuchos_ErrorPolling.hpp
  Teuchos_Exceptions.hpp
  Teuchos_F77_wrappers.h
  Teuchos_FancyOStream.hpp
  Teuchos_Flops.hpp
  Teuchos_GlobalMPISession.hpp
  Teuchos_Handle.hpp
  Teuchos_Handleable.hpp
  Teuchos_LAPACK.hpp
  Teuchos_LAPACK_wrappers.hpp
  Teuchos_LabeledObject.hpp
  Teuchos_MPIComm.hpp
  Teuchos_MPIContainerComm.hpp
  Teuchos_MPISession.hpp
  Teuchos_MPITraits.hpp
  Teuchos_NullIteratorTraits.hpp
  Teuchos_Object.hpp
  Teuchos_OpaqueWrapper.hpp
  Teuchos_OrdinalTraits.hpp
  Teuchos_ParameterEntry.hpp
  Teuchos_ParameterEntryValidator.hpp
  Teuchos_ParameterList.hpp
  Teuchos_ParameterListAcceptor.hpp
  Teuchos_ParameterListAcceptorDefaultBase.hpp
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
  Teuchos_SerialDenseMatrix.hpp
  Teuchos_SerialDenseVector.hpp
  Teuchos_SerializationTraits.hpp
  Teuchos_SerializationTraitsHelpers.hpp
  Teuchos_Serializer.hpp
  Teuchos_SerializerHelpers.hpp
  Teuchos_StandardCatchMacros.hpp
  Teuchos_StandardCompositionMacros.hpp
  Teuchos_StandardMemberCompositionMacros.hpp
  Teuchos_StandardParameterEntryValidators.hpp
  Teuchos_StringToIntMap.hpp
  Teuchos_toString.hpp
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
  Teuchos_UnitTestBase.hpp
  Teuchos_UnitTestHarness.hpp
  Teuchos_UnitTestHelpers.hpp
  Teuchos_UnitTestRepository.hpp
  Teuchos_Utils.hpp
  Teuchos_VerboseObject.hpp
  Teuchos_VerboseObjectParameterListHelpers.hpp
  Teuchos_VerbosityLevel.hpp
  Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp
  Teuchos_Workspace.hpp
  Teuchos_any.hpp
  Teuchos_arrayArg.hpp
  Teuchos_as.hpp
  Teuchos_basic_oblackholestream.hpp
  Teuchos_dyn_cast.hpp
  Teuchos_exit.h
  Teuchos_getConst.hpp
  Teuchos_implicit_cast.hpp
  Teuchos_map.hpp
  Teuchos_oblackholestream.hpp
  )

# Set sources

APPEND_SET(SOURCES
  Teuchos_BLAS.cpp
  Teuchos_CTimeMonitor.cpp
  Teuchos_CommandLineProcessor.cpp
  Teuchos_CompObject.cpp
  Teuchos_Describable.cpp
  Teuchos_ErrorPolling.cpp
  Teuchos_Flops.cpp
  Teuchos_GlobalMPISession.cpp
  Teuchos_LabeledObject.cpp
  Teuchos_MPIComm.cpp
  Teuchos_MPISession.cpp
  Teuchos_Object.cpp
  Teuchos_ParameterEntry.cpp
  Teuchos_ParameterList.cpp
  Teuchos_ParameterListAcceptor.cpp
  Teuchos_ParameterListAcceptorDefaultBase.cpp
  Teuchos_ParameterListNonAcceptor.cpp
  Teuchos_PerformanceMonitorUtils.cpp
  Teuchos_Ptr.cpp
  Teuchos_RCPNode.cpp
  Teuchos_Range1D.cpp
  Teuchos_ScalarTraits.cpp
  Teuchos_StandardParameterEntryValidators.cpp
  Teuchos_StringToIntMap.cpp
  Teuchos_TableColumn.cpp
  Teuchos_TableEntry.cpp
  Teuchos_TableFormat.cpp
  Teuchos_TabularOutputter.cpp
  Teuchos_TestForException.cpp
  Teuchos_Time.cpp
  Teuchos_TimeMonitor.cpp
  Teuchos_TypeNameTraits.cpp
  Teuchos_UnitTestBase.cpp
  Teuchos_UnitTestRepository.cpp
  Teuchos_Utils.cpp
  Teuchos_VerboseObject.cpp
  Teuchos_VerboseObjectParameterListHelpers.cpp
  Teuchos_VerbosityLevel.cpp
  Teuchos_VerbosityLevelCommandLineProcessorHelpers.cpp
  Teuchos_Workspace.cpp
  Teuchos_dyn_cast.cpp
  Teuchos_exit.cpp
  )

# Set addition headers and source based on package-specific
# configuration options

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)

  APPEND_SET(HEADERS
    Teuchos_ExpatHandlerAdapter.hpp
    Teuchos_FileInputSource.hpp
    Teuchos_FileInputStream.hpp
    Teuchos_HashSet.hpp
    Teuchos_Hashtable.hpp
    Teuchos_HashUtils.hpp
    Teuchos_Libxml2HandlerAdapter.hpp
    Teuchos_ParameterXMLFileReader.hpp
    Teuchos_StringInputSource.hpp
    Teuchos_StringInputStream.hpp
    Teuchos_StrUtils.hpp
    Teuchos_TreeBuildingXMLHandler.hpp
    Teuchos_XMLInputSource.hpp
    Teuchos_XMLInputStream.hpp
    Teuchos_XMLObject.hpp
    Teuchos_XMLObjectImplem.hpp
    Teuchos_XMLParameterListHelpers.hpp
    Teuchos_XMLParameterListReader.hpp
    Teuchos_XMLParameterListWriter.hpp
    Teuchos_XMLParser.hpp
    )

  APPEND_SET(SOURCES
    Teuchos_ExpatHandlerAdapter.cpp
    Teuchos_FileInputSource.cpp
    Teuchos_FileInputStream.cpp
    Teuchos_HashUtils.cpp
    Teuchos_Libxml2HandlerAdapter.cpp
    Teuchos_ParameterXMLFileReader.cpp
    Teuchos_StrUtils.cpp
    Teuchos_StringInputSource.cpp
    Teuchos_StringInputStream.cpp
    Teuchos_TreeBuildingXMLHandler.cpp
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

PACKAGE_ADD_LIBRARY(
  teuchos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
