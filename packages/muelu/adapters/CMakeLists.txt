SET(HEADERS "")
SET(SOURCES "")

#
# Package-specific configuration options
#

# access explicit instantiation
SET(BDIR ${CMAKE_CURRENT_BINARY_DIR})
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  INCLUDE_DIRECTORIES(${BDIR}/)
  APPEND_GLOB(HEADERS ${BDIR}/MueLu_ExplicitInstantiation.hpp)
ENDIF()

#
# Belos
#
IF (${PACKAGE_NAME}_ENABLE_Belos)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/belos)

  APPEND_SET(HEADERS
    belos/BelosOperatorT.hpp
    belos/BelosXpetraAdapter.hpp
    belos/BelosXpetraAdapterMultiVector.hpp
    belos/BelosXpetraAdapterOperator.hpp
    belos/BelosMueLuAdapter.hpp
    )

ENDIF()

#
# Epetra
#
IF (${PACKAGE_NAME}_ENABLE_Epetra)

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/epetra)

  APPEND_SET(HEADERS
    epetra/MueLu_EpetraOperator.hpp
    epetra/MueLu_CreateEpetraPreconditioner.hpp
    )

  APPEND_SET(SOURCES
    epetra/MueLu_EpetraOperator.cpp
    )

ENDIF()

#
# Tpetra
#
IF (${PACKAGE_NAME}_ENABLE_Tpetra)

   SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/tpetra)

   APPEND_SET(HEADERS
     tpetra/MueLu_CreateTpetraPreconditioner.hpp
     tpetra/MueLu_TpetraOperator_decl.hpp
     tpetra/MueLu_TpetraOperator_def.hpp
   )

   IF (${PACKAGE_NAME}_ENABLE_Belos AND TPETRA_INST_INT_INT)

     APPEND_SET(HEADERS
       tpetra/MueLu_RefMaxwell_decl.hpp
       tpetra/MueLu_RefMaxwell_def.hpp
       tpetra/MueLu_ShiftedLaplacian_decl.hpp
       tpetra/MueLu_ShiftedLaplacian_def.hpp
       tpetra/MueLu_ShiftedLaplacianOperator_decl.hpp
       tpetra/MueLu_ShiftedLaplacianOperator_def.hpp
     )

   ENDIF()

   TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

ENDIF()

#
# Stratimikos
#
IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stratimikos)
  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

  APPEND_SET(SOURCES
    ${DIR}/Stratimikos_MueLuHelpers.cpp
    )

  APPEND_SET(HEADERS
    ${DIR}/Stratimikos_MueLuHelpers.hpp
    ${DIR}/Thyra_MueLuPreconditionerFactory_decl.hpp
    ${DIR}/Thyra_MueLuPreconditionerFactory_def.hpp
    ${DIR}/Thyra_MueLuPreconditionerFactory.hpp
    )

  # Tpetra specific extensions
  IF (${PACKAGE_NAME}_ENABLE_Stratimikos AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
    SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/stratimikos)
    IF (${PACKAGE_NAME}_ENABLE_Teko AND ${PACKAGE_NAME}_ENABLE_Experimental)
      APPEND_SET(HEADERS
        ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory_decl.hpp
        ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory_def.hpp
        ${DIR}/Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp
        )
    ENDIF()
  ENDIF()
ENDIF()



IF (${PACKAGE_NAME}_ENABLE_AmgX)

  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/amgx)

  APPEND_GLOB(HEADERS ${DIR}/*.hpp)
  APPEND_GLOB(SOURCES ${DIR}/*.cpp)

  TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

ENDIF()

###################################################
## TAW: 9/22/15
## CMake-based ETI system inspired by Ifpack2 ETI
## - no support for KokkosClassic
## - ETI system for muelu-adapters only
## - The ETI system of MueLu is not touched!
###################################################

# Function to generate ETI (explicit template instantiation) files
# from a template and list of class names
FUNCTION(MUELUADAPTERS_PROCESS_ETI_TEMPLATE ETI_CLASSES TEMPLATE_FILE TPETRA_TEMPLATE_TYPE TEMPLATE_NODE PROCESSED_FILE SOURCES_LIST)
  SET(SRCS "")
  FOREACH(CLASS ${ETI_CLASSES})
    string(REPLACE "::" "_" CLASS_FILE_NAME "${CLASS}")
    string(TOUPPER "${CLASS_FILE_NAME}" UPPER_CASE_CLASS)
    string(REPLACE "CLASS_FILE_NAME" "${CLASS_FILE_NAME}" FINAL_FILE_NAME "${PROCESSED_FILE}")
    CONFIGURE_FILE(${TEMPLATE_FILE} ${FINAL_FILE_NAME})
    SET(SRCS ${SRCS} ${FINAL_FILE_NAME})
  ENDFOREACH()
  SET(${SOURCES_LIST} ${SRCS} PARENT_SCOPE)
ENDFUNCTION(MUELUADAPTERS_PROCESS_ETI_TEMPLATE)

SET(MUELUADAPTERS_ETI_CPP_SOURCES "")
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  # Set the list of classes in MueLu adapters that are templated
  # on <Scalar, LO, GO, Node> for which we want to do ETI using this system.
  IF (${PACKAGE_NAME}_ENABLE_Belos AND TPETRA_INST_INT_INT)
    GLOBAL_SET(MUELUADAPTERS_ETI_CLASSES
      MueLu::RefMaxwell
      )
  ENDIF()
  IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
    APPEND_GLOBAL_SET(MUELUADAPTERS_ETI_CLASSES
      Thyra::MueLuPreconditionerFactory
    )
  ENDIF()


  # Set the list of classes in MueLu adapters that are templated on <LO, GO, Node> for
  # which we want to do ETI using this system.
  GLOBAL_SET(MUELUADAPTERS_ETI_LO_GO_CLASSES
  )

  # ETI for the Kokkos::Serial device.
  IF (Tpetra_INST_SERIAL)

  ###  Kokkos_Compat_KokkosThreadsWrapperNode

    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_SC_LO_GO.tmpl
      "HAVE_TPETRA_INST_SERIAL"
      "Kokkos_Compat_KokkosSerialWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Serial.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_LO_GO_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_LO_GO.tmpl
      "HAVE_TPETRA_INST_SERIAL"
      "Kokkos_Compat_KokkosSerialWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Serial.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
  ENDIF()

  # ETI for the Kokkos::Threads device.
  IF (Tpetra_INST_PTHREAD)
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_SC_LO_GO.tmpl
      "HAVE_TPETRA_INST_PTHREAD"
      "Kokkos_Compat_KokkosThreadsWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Threads.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_LO_GO_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_LO_GO.tmpl
      "HAVE_TPETRA_INST_PTHREAD"
      "Kokkos_Compat_KokkosThreadsWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Threads.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
  ENDIF()

  # ETI for the Kokkos::OpenMP device.
  IF (Tpetra_INST_OPENMP)
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_SC_LO_GO.tmpl
      "HAVE_TPETRA_INST_OPENMP"
      "Kokkos_Compat_KokkosOpenMPWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_OpenMP.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_LO_GO_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_LO_GO.tmpl
      "HAVE_TPETRA_INST_OPENMP"
      "Kokkos_Compat_KokkosOpenMPWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_OpenMP.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
  ENDIF()

  # ETI for the Kokkos::Cuda device.
  IF (Tpetra_INST_CUDA)
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_SC_LO_GO.tmpl
      "HAVE_TPETRA_INST_CUDA"
      "Kokkos_Compat_KokkosCudaWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Cuda.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
    MUELUADAPTERS_PROCESS_ETI_TEMPLATE(
      "${MUELUADAPTERS_ETI_LO_GO_CLASSES}"
      ExplicitInstantiation/MueLu_ETI_LO_GO.tmpl
      "HAVE_TPETRA_INST_CUDA"
      "Kokkos_Compat_KokkosCudaWrapperNode"
      "ExplicitInstantiation/CLASS_FILE_NAME_Cuda.cpp"
      MueLuAdapters_SRCS)
    LIST(APPEND MUELUADAPTERS_ETI_CPP_SOURCES ${MueLuAdapters_SRCS})
  ENDIF()

  # add only automatically generated ETI cpp files
  APPEND_SET(SOURCES ${MUELUADAPTERS_ETI_CPP_SOURCES})
ENDIF()

############################################

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

TRIBITS_ADD_LIBRARY(
  muelu-adapters
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS muelu muelu-interface
  ADDED_LIB_TARGET_NAME_OUT MUELU_ADAPTERS_LIBNAME
  )

# Set linker language explicitly for CUDA builds.
set_target_properties(
  ${MUELU_ADAPTERS_LIBNAME}
  PROPERTIES LINKER_LANGUAGE CXX
  )
