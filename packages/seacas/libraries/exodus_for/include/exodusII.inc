                                ! -*- mode: fortran -*-
!
! Copyright(C) 1999-2020 National Technology & Engineering Solutions
! of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
! NTESS, the U.S. Government retains certain rights in this software.
!
! See packages/seacas/LICENSE for details
!***************************************************************
!
! Name: exodusII.inc - Exodus II Application include file
!
!       Use this include file to define general Exodus II parameters
!
!*****************************************************************************
!
!       Define Exodus functions that return values
!
        integer excre,exopen,exlgmd, exinqi
        external excre, exopen, exlgmd, exinqi
!
!       Define access modes
!
        integer EXNOCL,EXCLOB,EXREAD,EXWRIT,EXWRT,EXNORM,EXLARG
        integer EXNET4, EXNOSH, EXSHAR

        integer EX_READ,EX_WRITE,EX_NOCLOBBER,EX_CLOBBER,EX_NORMAL_MODEL
        integer EX_LARGE_MODEL,EX_NETCDF4,EX_NOSHARE,EX_SHARE
        integer EX_NOCLASSIC,EX_MAPS_INT64_DB,EX_IDS_INT64_DB
        integer EX_BULK_INT64_DB,EX_ALL_INT64_DB
        integer EX_MAPS_INT64_API
        integer EX_IDS_INT64_API,EX_BULK_INT64_API,EX_ALL_INT64_API
        integer EX_DISKLESS, EX_MMAP, EX_64BIT_OFFSET, EX_64BIT_DATA

        integer EX_OPT_MAX_NAME_LENGTH
        integer EX_OPT_COMPRESSION_TYPE
        integer EX_OPT_COMPRESSION_LEVEL
        integer EX_OPT_COMPRESSION_SHUFFLE
        integer EX_OPT_INTEGER_SIZE_API
        integer EX_OPT_INTEGER_SIZE_DB

        parameter (EXWRIT =     1)   !       EXWRIT - File open for write
        parameter (EXWRT  =     1)   !       EXWRT  - File open for write (alternate to EXWRIT)
        parameter (EXREAD =     2)   !       EXREAD - File open for read

        parameter (EXNOCL =     4)   !       EXNOCL - File create with no overwrite (no clobber), default
        parameter (EXCLOB =     8)   !       EXCLOB - File create with overwrite (clobber)
        parameter (EXNORM =    16)   !       EXNORM - Disable large model modifications
        parameter (EXLARG =    32)   !       EXLARG - Large model modifications
        parameter (EXNET4 =    64)   !       EXNET4 - use the hdf5-based netcdf4 output
        parameter (EXNOSH =   128)   !       EXNOSH - Do not open netcdf file in "share" mode
        parameter (EXSHAR =   256)   !       EXSHAR - Do open netcdf file in "share" mode

! Some of these are repeated from above, but use the same name as the C/C++ versions...

        parameter (EX_READ          =      0) !ex_open(): open file for reading (default)
        parameter (EX_WRITE         =      1) ! ex_open(): open existing file for appending.

        parameter (EX_NOCLOBBER     =      4) ! Don't overwrite existing database, default
        parameter (EX_CLOBBER       =      8) ! Overwrite existing database if it exists
        parameter (EX_NORMAL_MODEL  =     16) ! disable mods that permit storage of larger models
        parameter (EX_LARGE_MODEL   =     32) ! enable mods that permit storage of larger models
        parameter (EX_NETCDF4       =     64) ! use the hdf5-based netcdf4 output
        parameter (EX_64BIT_OFFSET  =     32) ! enable mods that permit storage of larger models
        parameter (EX_64BIT_DATA    =4194304) ! CDF-5 format: classic model but 64 bit dimensions and sizes
        parameter (EX_NOSHARE       =    128) ! Do not open netcdf file in "share" mode
        parameter (EX_SHARE         =    256) ! Do open netcdf file in "share" mode
        parameter (EX_NOCLASSIC     =    512) ! Do not force netcdf to classic mode in netcdf4 mode
        parameter (EX_DISKLESS      =1048576) ! Experimental
        parameter (EX_MMAP          =2097152) ! Experimental

! Need to distinguish between storage on database (DB in name) and
! passed through the API functions (API in name).

        parameter (EX_MAPS_INT64_DB =   1024) ! All maps (id, order, ...) store int64_t values
        parameter (EX_IDS_INT64_DB  =   2048) ! All entity ids (sets, blocks, maps) are int64_t values
        parameter (EX_BULK_INT64_DB =   4096) ! All integer bulk data (local indices, counts, maps); not ids
        parameter (EX_ALL_INT64_DB  =   1024+2048+4096) ! All of the above...

        parameter (EX_MAPS_INT64_API=   8192) ! All maps (id, order, ...) store int64_t values
        parameter (EX_IDS_INT64_API =  16384) ! All entity ids (sets, blocks, maps) are int64_t values
        parameter (EX_BULK_INT64_API=  32768) ! All integer bulk data (local indices, counts, maps); not ids
        parameter (EX_ALL_INT64_API =  8192+16384+32768) ! All of the above...

!
!  -- For EXSETOPT Calls
        parameter (EX_OPT_MAX_NAME_LENGTH = 1)
        parameter (EX_OPT_COMPRESSION_TYPE = 2)    !Currently not used. GZip by default
        parameter (EX_OPT_COMPRESSION_LEVEL = 3)   ! 0 (disabled/fastest) ... 9 (best/slowest)
        parameter (EX_OPT_COMPRESSION_SHUFFLE = 4) ! 0 (disabled); 1 (enabled)
        parameter (EX_OPT_INTEGER_SIZE_API = 5)    ! See *_INT64_* values above
        parameter (EX_OPT_INTEGER_SIZE_DB  = 6)    ! (query only)
!
!  --INQUIRE FLAGS
        integer EX_INQ_FILE_TYPE, EX_INQ_API_VERS, EX_INQ_DB_VERS
        integer EX_INQ_TITLE, EX_INQ_DIM, EX_INQ_NODES
        integer EX_INQ_ELEM, EX_INQ_ELEM_BLK, EX_INQ_NODE_SETS
        integer EX_INQ_NS_NODE_LEN, EX_INQ_SIDE_SETS
        integer EX_INQ_SS_NODE_LEN, EX_INQ_SS_ELEM_LEN, EX_INQ_QA
        integer EX_INQ_INFO, EX_INQ_TIME, EX_INQ_EB_PROP
        integer EX_INQ_NS_PROP, EX_INQ_SS_PROP, EX_INQ_NS_DF_LEN
        integer EX_INQ_SS_DF_LEN, EX_INQ_LIB_VERS
        integer EX_INQ_EM_PROP, EX_INQ_NM_PROP, EX_INQ_ELEM_MAP
        integer EX_INQ_NODE_MAP, EX_INQ_EDGE
        integer EX_INQ_EDGE_BLK, EX_INQ_EDGE_SETS, EX_INQ_ES_LEN
        integer EX_INQ_ES_DF_LEN, EX_INQ_EDGE_PROP
        integer EX_INQ_ES_PROP, EX_INQ_FACE, EX_INQ_FACE_BLK
        integer EX_INQ_FACE_SETS, EX_INQ_FS_LEN
        integer EX_INQ_FS_DF_LEN, EX_INQ_FACE_PROP, EX_INQ_FS_PROP
        integer EX_INQ_ELEM_SETS, EX_INQ_ELS_LEN
        integer EX_INQ_ELS_DF_LEN, EX_INQ_ELS_PROP, EX_INQ_EDGE_MAP
        integer EX_INQ_FACE_MAP, EX_INQ_COORD_FRAMES
        integer EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH
        integer EX_INQ_DB_MAX_USED_NAME_LENGTH
        integer EX_INQ_MAX_READ_NAME_LENGTH
        integer EX_INQ_DB_FLOAT_SIZE, EX_INQ_NUM_CHILD_GROUPS
        integer EX_INQ_GROUP_PARENT, EX_INQ_GROUP_ROOT
        integer EX_INQ_GROUP_NAME_LEN, EX_INQ_GROUP_NAME
        integer EX_INQ_FULL_GROUP_NAME_LEN, EX_INQ_FULL_GROUP_NAME
        integer EX_INQ_THREADSAFE, EX_INQ_INVALID

        parameter (EX_INQ_FILE_TYPE    = 1) !   inquire EXODUS file type
        parameter (EX_INQ_API_VERS     = 2) !   inquire API version number
        parameter (EX_INQ_DB_VERS      = 3) !   inquire database version number
        parameter (EX_INQ_TITLE        = 4) !   inquire database title
        parameter (EX_INQ_DIM          = 5) !   inquire number of dimensions
        parameter (EX_INQ_NODES        = 6) !   inquire number of nodes
        parameter (EX_INQ_ELEM         = 7) !   inquire number of elements
        parameter (EX_INQ_ELEM_BLK     = 8) !   inquire number of element blocks
        parameter (EX_INQ_NODE_SETS    = 9) !   inquire number of node sets
        parameter (EX_INQ_NS_NODE_LEN  = 10) !  inquire length of node set node list
        parameter (EX_INQ_SIDE_SETS    = 11) !  inquire number of side sets
        parameter (EX_INQ_SS_NODE_LEN  = 12) !  inquire length of side set node list
        parameter (EX_INQ_SS_ELEM_LEN  = 13) !  inquire length of side set element list
        parameter (EX_INQ_QA           = 14) !  inquire number of QA records
        parameter (EX_INQ_INFO         = 15) !  inquire number of info records
        parameter (EX_INQ_TIME         = 16) !  inquire number of time steps in the database
        parameter (EX_INQ_EB_PROP      = 17) !  inquire number of element block properties
        parameter (EX_INQ_NS_PROP      = 18) !  inquire number of node set properties
        parameter (EX_INQ_SS_PROP      = 19) !  inquire number of side set properties
        parameter (EX_INQ_NS_DF_LEN    = 20) !  inquire length of node set distribution factor list
        parameter (EX_INQ_SS_DF_LEN    = 21) !  inquire length of side set distribution factor list
        parameter (EX_INQ_LIB_VERS     = 22) !  inquire API Lib vers number
        parameter (EX_INQ_EM_PROP      = 23) !  inquire number of element map properties
        parameter (EX_INQ_NM_PROP      = 24) !  inquire number of node map properties
        parameter (EX_INQ_ELEM_MAP     = 25) !  inquire number of element maps
        parameter (EX_INQ_NODE_MAP     = 26) !  inquire number of node maps
        parameter (EX_INQ_EDGE         = 27) !  inquire number of edges
        parameter (EX_INQ_EDGE_BLK     = 28) !  inquire number of edge blocks
        parameter (EX_INQ_EDGE_SETS    = 29) !  inquire number of edge sets
        parameter (EX_INQ_ES_LEN       = 30) !  inquire length of concat edge set edge list
        parameter (EX_INQ_ES_DF_LEN    = 31) !  inquire length of concat edge set dist factor list
        parameter (EX_INQ_EDGE_PROP    = 32) !  inquire number of properties stored per edge block
        parameter (EX_INQ_ES_PROP      = 33) !  inquire number of properties stored per edge set
        parameter (EX_INQ_FACE         = 34) !  inquire number of faces
        parameter (EX_INQ_FACE_BLK     = 35) !  inquire number of face blocks
        parameter (EX_INQ_FACE_SETS    = 36) !  inquire number of face sets
        parameter (EX_INQ_FS_LEN       = 37) !  inquire length of concat face set face list
        parameter (EX_INQ_FS_DF_LEN    = 38) !  inquire length of concat face set dist factor list
        parameter (EX_INQ_FACE_PROP    = 39) !  inquire number of properties stored per face block
        parameter (EX_INQ_FS_PROP      = 40) !  inquire number of properties stored per face set
        parameter (EX_INQ_ELEM_SETS    = 41) !  inquire number of element sets
        parameter (EX_INQ_ELS_LEN      = 42) !  inquire length of concat element set element list
        parameter (EX_INQ_ELS_DF_LEN   = 43) !  inquire length of concat element set dist factor list
        parameter (EX_INQ_ELS_PROP     = 44) !  inquire number of properties stored per elem set
        parameter (EX_INQ_EDGE_MAP     = 45) !  inquire number of edge maps
        parameter (EX_INQ_FACE_MAP     = 46) !  inquire number of face maps
        parameter (EX_INQ_COORD_FRAMES = 47) !  inquire number of coordinate frames
        parameter (EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH = 48) !inquire size of MAX_NAME_LENGTH dimension on database
        parameter (EX_INQ_DB_MAX_USED_NAME_LENGTH = 49) !  inquire size of MAX_NAME_LENGTH dimension on database
        parameter (EX_INQ_MAX_READ_NAME_LENGTH    = 50) !  inquire client-specified max size of returned names
        parameter (EX_INQ_DB_FLOAT_SIZE           = 51) !  inquire size of floating-point values stored on database
        parameter (EX_INQ_NUM_CHILD_GROUPS        = 52) !  inquire number of groups contained in this (exoid) group
        parameter (EX_INQ_GROUP_PARENT            = 53) !  inquire id of parent of this (exoid) group; returns exoid if at root
        parameter (EX_INQ_GROUP_ROOT              = 54) !  inquire id of root group "/" of this (exoid) group; returns exoid if at root
        parameter (EX_INQ_GROUP_NAME_LEN          = 55) !  inquire length of name of group exoid
        parameter (EX_INQ_GROUP_NAME              = 56) !  inquire name of group exoid. "/" returned for root group
        parameter (EX_INQ_FULL_GROUP_NAME_LEN     = 57) !  inquire length of full path name of this (exoid) group
        parameter (EX_INQ_FULL_GROUP_NAME         = 58) !  inquire full "/"-separated path name of this (exoid) group
        parameter (EX_INQ_THREADSAFE              = 59) !  Returns 1 if library is thread-safe; 0 otherwise
        parameter (EX_INQ_INVALID                 = -1)

      integer EXTYPE, EXVERS, EXDBVR, EXTITL, EXDIM, EXNODE
      integer  EXELEM, EXELBL, EXNODS, EXSIDS, EXQA,  EXINFO
      integer  EXTIMS, EXSSNL, EXSSEL, EXNSNL, EXNEBP, EXNNSP
      integer  EXNSSP, EXNSDF, EXSSDF, EXLBVR, EXNEMP, EXNNMP
      integer  EXNEM,  EXNNM,  EXNEDG, EXNEDB, EXNEDS, EXESL
      integer  EXESDL, EXEDBP, EXEDSP, EXNFAC, EXNFB , EXNFS
      integer  EXFSL , EXFSDL, EXFABP, EXFASP, EXNELS, EXELSL
      integer  EXELSD, EXELSP, EXNEDM, EXNFAM, EXNCF,  EXDBMXALNM
      integer  EXDBMXUSNM, EXMXRDNM

!
      parameter (EXTYPE = 1)     ! --  inquire EXODUS II file type
      parameter (EXVERS = 2)     ! --  inquire API version number
      parameter (EXDBVR = 3)     ! --  inquire database version number
      parameter (EXTITL = 4)     ! --  inquire database title
      parameter (EXDIM  = 5)     ! --  inquire number of dimensions
      parameter (EXNODE = 6)     ! --  inquire number of nodes
      parameter (EXELEM = 7)     ! --  inquire number of elements
      parameter (EXELBL = 8)     ! --  inquire number of element blocks
      parameter (EXNODS = 9)     ! --  inquire number of node sets
      parameter (EXNSNL = 10)    ! --  inquire length of node set node list
      parameter (EXSIDS = 11)    ! --  inquire number of side sets
      parameter (EXSSNL = 12)    ! --  inquire length of side set node list
      parameter (EXSSEL = 13)    ! --  inquire length of side set element list
      parameter (EXQA   = 14)    ! --  inquire number of QA records
      parameter (EXINFO = 15)    ! --  inquire number of info records
      parameter (EXTIMS = 16)    ! --  inquire number of time steps in the database
      parameter (EXNEBP = 17)    ! --  inquire number of element block properties
      parameter (EXNNSP = 18)    ! --  inquire number of node set properties
      parameter (EXNSSP = 19)    ! --  inquire number of side set properties
      parameter (EXNSDF = 20)    ! --  inquire length of node set distribution factor list
      parameter (EXSSDF = 21)    ! --  inquire length of side set distribution factor list
      parameter (EXLBVR = 22)    ! --  inquire API Lib vers number
      parameter (EXNEMP = 23)    ! --  inquire number of element map properties
      parameter (EXNNMP = 24)    ! --  inquire number of node map properties
      parameter (EXNEM  = 25)    ! --  inquire number of element maps
      parameter (EXNNM  = 26)    ! --  inquire number of node maps
      parameter (EXNEDG = 27)    ! --  inquire number of edges
      parameter (EXNEDB = 28)    ! --  inquire number of edge blocks
      parameter (EXNEDS = 29)    ! --  inquire number of edge sets
      parameter (EXESL  = 30)    ! --  inquire length of concat edge set edge list
      parameter (EXESDL = 31)    ! --  inquire length of concat edge set dist factor list
      parameter (EXEDBP = 32)    ! --  inquire number of properties stored per edge block
      parameter (EXEDSP = 33)    ! --  inquire number of properties stored per edge set
      Parameter (EXNFAC = 34)    ! --  inquire number of faces
      parameter (EXNFB  = 35)    ! --  inquire number of face blocks
      parameter (EXNFS  = 36)    ! --  inquire number of face sets
      parameter (EXFSL  = 37)    ! --  inquire length of concat face set face list
      parameter (EXFSDL = 38)    ! --  inquire length of concat face set dist factor list
      parameter (EXFABP = 39)    ! --  inquire number of properties stored per face block
      parameter (EXFASP = 40)    ! --  inquire number of properties stored per face set
      parameter (EXNELS = 41)    ! --  inquire number of element sets
      parameter (EXELSL = 42)    ! --  inquire length of concat element set element list
      parameter (EXELSD = 43)    ! --  inquire length of concat element set dist factor list
      parameter (EXELSP = 44)    ! --  inquire number of properties stored per elem set
      parameter (EXNEDM = 45)    ! --  inquire number of edge maps
      parameter (EXNFAM = 46)    ! --  inquire number of face maps
      parameter (EXNCF  = 47)    ! --  inquire number of coordinate frames
      parameter (EXDBMXALNM = 48)! --  inquire size of MAX_NAME_LENGTH dimension on database
      parameter (EXDBMXUSNM = 49)! --  inquire max size of any name existing on database
      parameter (EXMXRDNM = 50)  ! --  inquire user-specified max size of names returned; truncated if longer
!
!         Define limits
!
      integer MXSTLN, MXLNLN, MXNAME
!     MXSTLN - Maximum string length: used for store strings
      parameter (MXSTLN = 32)
!     MXLNLN - Maximum line string length
      parameter (MXLNLN = 80)
!     MXNAME - Maximum possible length of a name (due to netcdf NC_MAX_NAME)
      parameter (MXNAME = 256)
!
!     properties
!
      integer EXEBLK, EXNSET, EXSSET, EXEMAP, EXNMAP

      parameter (EXEBLK = 1)   ! EXEBLK - element block property code
      parameter (EXNSET = 2)   ! EXNSPR - node set property code
      parameter (EXSSET = 3)   ! EXSSPR - side set property code
      parameter (EXEMAP = 4)   ! EXEMAP - element map property code
      parameter (EXNMAP = 5)   ! EXNMAP - node map property code

!     Define values for type
!

      integer EX_NODAL, EX_NODE_BLOCK, EX_NODE_SET, EX_EDGE_BLOCK
      integer EX_EDGE_SET, EX_FACE_BLOCK, EX_FACE_SET, EX_ELEM_BLOCK
      integer EX_ELEM_SET, EX_SIDE_SET, EX_ELEM_MAP, EX_NODE_MAP
      integer EX_EDGE_MAP, EX_FACE_MAP, EX_GLOBAL, EX_COORDINATE
      integer EX_INVALID

      parameter (EX_NODAL       = 14)      !     nodal "block" for variables
      parameter (EX_NODE_BLOCK  = 14)      !     alias for EX_NODAL
      parameter (EX_NODE_SET    =  2)      !     node set property code
      parameter (EX_EDGE_BLOCK  =  6)      !     edge block property code
      parameter (EX_EDGE_SET    =  7)      !     edge set property code
      parameter (EX_FACE_BLOCK  =  8)      !     face block property code
      parameter (EX_FACE_SET    =  9)      !     face set property code
      parameter (EX_ELEM_BLOCK  =  1)      !     element block property code
      parameter (EX_ELEM_SET    = 10)      !     face set property code
      parameter (EX_SIDE_SET    =  3)      !     side set property code

      parameter (EX_ELEM_MAP    =  4)      !     element map property code
      parameter (EX_NODE_MAP    =  5)      !     node map property code
      parameter (EX_EDGE_MAP    = 11)      !     edge map property code
      parameter (EX_FACE_MAP    = 12)      !     face map property code

      parameter (EX_GLOBAL      = 13)      !     global "block" for variables
      parameter (EX_COORDINATE  = 15)      !     kluge so some internal wrapper functions work
      parameter (EX_INVALID     = -1)
!
!     Define options codes
!
      integer EXVRBS, EXDEBG, EXABRT

      parameter (EXVRBS = 1)   ! EXVRBS - verbose mode message flag
      parameter (EXDEBG = 2)   ! EXDEBG - debug mode message flag
      parameter (EXABRT = 4)   ! EXABRT - abort mode flag
!
!     Define coordinate frames
!
      integer EXCFREC, EXCFCYL, EXCFSPH
      parameter (EXCFREC = 1)     ! Rectangular
      parameter (EXCFCYL = 2)     ! Cylindrical
      parameter (EXCFSPH = 3)     ! Spherical

!
!     Define error return codes
!
      integer EXFATL, EXOK, EXWARN

      parameter (EXFATL = -1)  ! EXFATL - fatal error flag def
      parameter (EXOK = 0)     ! EXOK - no err (success) flag def
      parameter (EXWARN = 1)   ! EXWARN - warning flag def
!
!     Define error return codes and flags
!
      integer EXMEMF, EXBFMD, EXBFID, EXBTID, EXBPRM, EXNULL
      integer EXPMSG, EXLMSG

      parameter (EXMEMF =  1000)   ! EXMEMF - memory allocation failure flag def
      parameter (EXBFMD =  1001)   ! EXBFMD - wrong file mode
      parameter (EXBFID =  1002)   ! EXBFID - bad file id def
      parameter (EXBTID =  1004)   ! EXBTID - id table lookup failed
      parameter (EXBPRM =  1005)   ! EXBPRM - bad parameter passed
      parameter (EXNULL = -1006)   ! EXNULL - null entity found
      parameter (EXPMSG = -1000)   ! EXPMSG - message print code - no error implied
      parameter (EXLMSG = -1001)   ! EXLMSG - print last error message msg code
