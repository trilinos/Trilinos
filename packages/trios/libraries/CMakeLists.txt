
INCLUDE(PackageMacros)

ADD_SUBDIRECTORY(support)
#ADD_SUBDIRECTORY(commsplitter)

# We need to make sure at least one of the optional TPLs for data transport
# was enabled. There might be a way to do this in the Dependencies.cmake file, 
# but I couldn't figure it out. 

IF (${PACKAGE_NAME}_ENABLE_Portals OR
    ${PACKAGE_NAME}_ENABLE_CrayPortals OR
    ${PACKAGE_NAME}_ENABLE_InfiniBand OR
    ${PACKAGE_NAME}_ENABLE_Gemini OR 
    ${PACKAGE_NAME}_ENABLE_LUC)
        
    ADD_SUBDIRECTORY(nessie)
    
    # If Gemini is the transport, we need the commsplitter lib as well
    IF (${PACKAGE_NAME}_ENABLE_Gemini)
      ADD_SUBDIRECTORY(commsplitter)
    ENDIF()
    
ELSE ()
    MESSAGE(WARNING "Did not find a supported network protocol. "
            "Nessie requires one of the following TPLs to be enabled: "
            " Portals, CrayPortals, InfiniBand, LUC, Gemini")
ENDIF()

