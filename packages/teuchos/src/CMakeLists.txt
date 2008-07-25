# $Header$ 

# Package-specific configuration options
IF(${PROJECT_NAME}_ENABLE_COMPLEX)
  SET(HAVE_TEUCHOS_COMPLEX TRUE)
ENDIF(${PROJECT_NAME}_ENABLE_COMPLEX)

IF(${PROJECT_NAME}_ENABLE_EXTENDED)
  SET(HAVE_TEUCHOS_EXTENDED TRUE)
ENDIF(${PROJECT_NAME}_ENABLE_EXTENDED)

# Create <PACKAGE>_config.h
CONFIGURE_FILE(
  ${${PROJECT_NAME}_SOURCE_DIR}/cmake/${PROJECT_NAME}_config.h.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h)
INSTALL(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_config.h
  DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Setup directories based on current directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Set headers
SET(HEADERS
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
	#Teuchos_PrimitiveTypeTraits.hpp 
	Teuchos_PromotionTraits.hpp 
	Teuchos_Ptr.hpp 
	Teuchos_RCP.hpp 
	Teuchos_RCPNode.hpp 
	Teuchos_RCPBoostSharedPtrConversions.hpp 
	Teuchos_RCPBoostSharedPtrConversionsDecl.hpp 
	Teuchos_RCPDecl.hpp 
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
	Teuchos_TableColumn.hpp 
	Teuchos_TableEntry.hpp 
	Teuchos_TableFormat.hpp 
	Teuchos_TestForException.hpp 
	Teuchos_TestingHelpers.hpp 
	Teuchos_Time.hpp 
	Teuchos_TimeMonitor.hpp 
	Teuchos_TypeNameTraits.hpp 
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
SET(SOURCES
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
	Teuchos_TestForException.cpp 
	Teuchos_Time.cpp 
	Teuchos_TimeMonitor.cpp 
	Teuchos_TypeNameTraits.cpp 
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
IF(${PROJECT_NAME}_ENABLE_EXTENDED)
	SET(HEADERS
		${HEADERS}
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

	SET(SOURCES
		${SOURCES}
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

ENDIF(${PROJECT_NAME}_ENABLE_EXTENDED)

IF(TRILINOS_ENABLE_MPI)

  SET(HEADERS
    Teuchos_DefaultMpiComm.hpp
    Teuchos_MpiReductionOpSetter.hpp
    ${HEADERS}
    )

  SET(SOURCES
    ${SOURCES}
    Teuchos_MpiReductionOpSetter.cpp
    )
  
ENDIF(TRILINOS_ENABLE_MPI)

# 2008/06/09: rabartl: Above: We need to change all of this to use file
# (GLOB...) and then use ifdefs in the files to control what gets built.

# 2008/06/11: dmdunla: file globbing will work if the directories contain
# only those files that should be built. Currently, though, there are 
# extraneous files (e.g., Teuchos_FloatingPointTrap.cpp)

# Target library setup
ADD_LIBRARY(teuchos ${HEADERS} ${SOURCES})
LIST(APPEND TARGET_LIBS teuchos)
LIST(APPEND LIB_DEPS ${TRILINOS_LAPACK_LIBRARY} ${TRILINOS_BLAS_LIBRARY})
TARGET_LINK_LIBRARIES(teuchos ${LIB_DEPS})

# 2008/06/09: rabartl: Above: Why are the headers being listed in
# ADD_LIBRARY(...)?  Can't Cmake figure out what the headers are?

# 2008/06/11: dmdunla: The addition of headers in the targets
# enables access to headers from within MS VisualStudio

# MPI setup
IF(TRILINOS_ENABLE_MPI)
  LIST(APPEND LIB_DEPS ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(teuchos ${LIB_DEPS})
ENDIF(TRILINOS_ENABLE_MPI)

INSTALL(
  TARGETS teuchos
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib
  COMPONENT ${PROJECT_NAME}
	)

INSTALL(
  FILES ${HEADERS}
  DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR}
  COMPONENT ${PROJECT_NAME}
  )

# Define current package dependency variables
GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
SET(${PROJECT_NAME}_INCLUDE_DIRS ${${PROJECT_NAME}_INCLUDE_DIRS} ${INCLUDE_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} include directories")
GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)
SET(${PROJECT_NAME}_LIBRARY_DIRS ${${PROJECT_NAME}_LIBRARY_DIRS}  ${LIBRARY_DIRS_CURRENT} CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARIES} ${TARGET_LIBS} ${LIB_DEPS} CACHE INTERNAL "${PROJECT_NAME} libraries")
