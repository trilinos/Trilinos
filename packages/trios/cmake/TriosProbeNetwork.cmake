
include(CheckLibraryExists)
include(CheckFunctionExists)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckTypeSize)


########## INCLUDE FILES ##############
include(CheckIncludeFiles)

check_include_files("malloc.h" HAVE_TRIOS_MALLOC_H)
check_include_files("sys/types.h" HAVE_TRIOS_SYS_TYPES_H)
check_include_files("sys/param.h" HAVE_TRIOS_SYS_PARAM_H)
check_include_files("sys/ioctl.h" HAVE_TRIOS_SYS_IOCTL_H)
check_include_files("sys/socket.h" HAVE_TRIOS_SYS_SOCKET_H)
check_include_files("sys/sockio.h" HAVE_TRIOS_SYS_SOCKIO_H)
check_include_files("netdb.h" HAVE_TRIOS_NETDB_H)
check_include_files("sys/socket;net/if.h" HAVE_TRIOS_NET_IF_H)
check_include_files("sys/socket.h;net/if_dl.h" HAVE_TRIOS_NET_IF_DL_H)
check_include_files("sys/socket.h;net/if_arp.h" HAVE_TRIOS_NET_IF_ARP_H)
check_include_files("netinet/in.h" HAVE_TRIOS_NETINET_IN_H)
check_include_files("arpa/inet.h" HAVE_TRIOS_ARPA_INET_H)
check_include_files("sys/types.h;sys/socket.h;ifaddrs.h" HAVE_TRIOS_IFADDRS_H)

########## Probe for various network configurations ##############


check_function_exists(getifaddrs TRIOS_HAVE_GETIFADDRS)

set(CMAKE_EXTRA_INCLUDE_FILES "netinet/in.h")
check_type_size("struct sockaddr_in" HAVE_TRIOS_STRUCT_SOCKADDR_IN)
set(CMAKE_EXTRA_INCLUDE_FILES)


######## Portals Config ##################

IF (${PACKAGE_NAME}_ENABLE_Portals OR ${PACKAGE_NAME}_ENABLE_CrayPortals)

    message(STATUS "Checking Portals configuration: PORTALS_IFACE_DEFAULT=${PORTALS_IFACE_DEFAULT}")

    # Need this macro defined to get the right integer types for Portals
    ADD_DEFINITIONS(-D__STDC_CONSTANT_MACROS)

    # Also need values for PTL_IFACE_SERVER and PTL_IFACE_CLIENT
    if (NOT PORTALS_IFACE_SERVER)
        set(PORTALS_IFACE_SERVER "PTL_IFACE_DEFAULT" CACHE STRING "Default NAL for the Portals server")
    endif (NOT PORTALS_IFACE_SERVER)
        if (NOT PORTALS_IFACE_CLIENT)
        set(PORTALS_IFACE_CLIENT "PTL_IFACE_DEFAULT" CACHE STRING "Default NAL for a Portals client")
    endif (NOT PORTALS_IFACE_CLIENT)

    ADD_DEFINITIONS(-DPTL_IFACE_CLIENT=${PORTALS_IFACE_CLIENT})
    ADD_DEFINITIONS(-DPTL_IFACE_SERVER=${PORTALS_IFACE_SERVER})

    # Figure out the name of the Portals include file
    IF (${PACKAGE_NAME}_ENABLE_Portals)
      SET (PORTALS_INCLUDE_FILE "portals3.h")
    ENDIF (${PACKAGE_NAME}_ENABLE_Portals)

    IF (${PACKAGE_NAME}_ENABLE_CrayPortals)
      SET (PORTALS_INCLUDE_FILE "portals/portals3.h")
    ENDIF (${PACKAGE_NAME}_ENABLE_CrayPortals)

    SET(CMAKE_REQUIRED_INCLUDES ${Portals_INCLUDE_DIRS})
    SET(CMAKE_EXTRA_INCLUDE_FILES ${PORTALS_INCLUDE_FILE})
    set(CMAKE_REQUIRED_LIBRARIES ${TPL_Portals_LIBRARIES})

    check_include_files(p3nal_utcp.h HAVE_TRIOS_P3NAL_UTCP_H)
    check_include_files("p3nal_utcp.h" HAVE_TRIOS_P3NAL_UTCP_H)
    check_include_files("p3rt/p3rt.h" HAVE_TRIOS_P3RT_P3RT_H)

    # Check for PTL_NOACK_REQ
    check_c_source_compiles(
        "#include<${PORTALS_INCLUDE_FILE}>\nmain(){return PTL_NOACK_REQ;}"
        HAVE_TRIOS_PTL_NOACK_REQ)

    # Check for PTL_NO_ACK_REQ
    check_c_source_compiles(
        "#include<${PORTALS_INCLUDE_FILE}>\nint main() {return PTL_NO_ACK_REQ;}"
        HAVE_TRIOS_PTL_NO_ACK_REQ)

    # Check for portals types

    check_type_size(ptl_time_t HAVE_TRIOS_PTL_TIME_T)
    check_type_size(ptl_eq_handler_t HAVE_TRIOS_PTL_EQ_HANDLER_T)

    # Check for portals functions
    check_function_exists(PtlErrorStr HAVE_TRIOS_PTLERRORSTR)
    check_function_exists(PtlNIFailStr HAVE_TRIOS_PTLNIFAILSTR)
    check_function_exists(PtlEventKindStr HAVE_TRIOS_PTLEVENTKINDSTR)
    check_function_exists(PtlGetJid HAVE_TRIOS_PTLGITJID)
    check_function_exists(PtlACEntry HAVE_TRIOS_PTLACENTRY)

    SET(CMAKE_REQUIRED_INCLUDES)
    SET(CMAKE_EXTRA_INCLUDE_FILES)
    set(CMAKE_REQUIRED_LIBRARIES)

ENDIF (${PACKAGE_NAME}_ENABLE_Portals OR ${PACKAGE_NAME}_ENABLE_CrayPortals)


