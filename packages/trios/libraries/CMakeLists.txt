
INCLUDE(TribitsPackageMacros)

# We need to make sure at least one of the optional TPLs for data transport
# was enabled. There might be a way to do this in the Dependencies.cmake file,
# but I couldn't figure it out.

IF (${PACKAGE_NAME}_HAVE_RDMA_TRANSPORT)

  # If Gemini is the transport, we need the commsplitter lib so we
  # can run services in MPMD mode
  IF (${PACKAGE_NAME}_ENABLE_Gemini)
    ADD_SUBDIRECTORY(commsplitter)
  ENDIF()

ELSE ()
    MESSAGE(WARNING "Did not find a supported network protocol. "
            "Nessie requires one of the following TPLs to be enabled: "
            " Portals, CrayPortals, InfiniBand, LUC, Gemini")
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/nessie/client)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/nessie/server)
