! -*- mode: fortran -*-
!
! Copyright(C) 1999-2020, 2025 National Technology & Engineering Solutions
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
        integer, parameter :: EXWRIT =     1 !       EXWRIT - File open for write
        integer, parameter :: EXWRT  =     1 !       EXWRT  - File open for write (alternate to EXWRIT)
        integer, parameter :: EXREAD =     2 !       EXREAD - File open for read

        integer, parameter :: EXNOCL =     4 !       EXNOCL - File create with no overwrite (no clobber), default
        integer, parameter :: EXCLOB =     8 !       EXCLOB - File create with overwrite (clobber)
        integer, parameter :: EXNORM =    16 !       EXNORM - Disable large model modifications
        integer, parameter :: EXLARG =    32 !       EXLARG - Large model modifications
        integer, parameter :: EXNET4 =    64 !       EXNET4 - use the hdf5-based netcdf4 output
        integer, parameter :: EXNOSH =   128 !       EXNOSH - Do not open netcdf file in "share" mode
        integer, parameter :: EXSHAR =   256 !       EXSHAR - Do open netcdf file in "share" mode

! Some of these are repeated from above, but use the same name as the C/C++ versions...
        integer, parameter :: EX_READ          =      0 !ex_open(): open file for reading (default)
        integer, parameter :: EX_WRITE         =      1 ! ex_open(): open existing file for appending.

        integer, parameter :: EX_NOCLOBBER     =      4 ! Don't overwrite existing database, default
        integer, parameter :: EX_CLOBBER       =      8 ! Overwrite existing database if it exists
        integer, parameter :: EX_NORMAL_MODEL  =     16 ! disable mods that permit storage of larger models
        integer, parameter :: EX_LARGE_MODEL   =     32 ! enable mods that permit storage of larger models
        integer, parameter :: EX_NETCDF4       =     64 ! use the hdf5-based netcdf4 output
        integer, parameter :: EX_64BIT_OFFSET  =     32 ! enable mods that permit storage of larger models
        integer, parameter :: EX_64BIT_DATA    =4194304 ! CDF-5 format: classic model but 64 bit dimensions and sizes
        integer, parameter :: EX_NOSHARE       =    128 ! Do not open netcdf file in "share" mode
        integer, parameter :: EX_SHARE         =    256 ! Do open netcdf file in "share" mode
        integer, parameter :: EX_NOCLASSIC     =    512 ! Do not force netcdf to classic mode in netcdf4 mode
        integer, parameter :: EX_DISKLESS      =1048576 ! Experimental
        integer, parameter :: EX_MMAP          =2097152 ! Experimental

! Need to distinguish between storage on database (DB in name) and
! passed through the API functions (API in name).
        integer, parameter :: EX_MAPS_INT64_DB =   1024 ! All maps (id, order, ...) store int64_t values
        integer, parameter :: EX_IDS_INT64_DB  =   2048 ! All entity ids (sets, blocks, maps) are int64_t values
        integer, parameter :: EX_BULK_INT64_DB =   4096 ! All integer bulk data (local indices, counts, maps); not ids
        integer, parameter :: EX_ALL_INT64_DB  =   1024+2048+4096 ! All of the above...

        integer, parameter :: EX_MAPS_INT64_API=   8192 ! All maps (id, order, ...) store int64_t values
        integer, parameter :: EX_IDS_INT64_API =  16384 ! All entity ids (sets, blocks, maps) are int64_t values
        integer, parameter :: EX_BULK_INT64_API=  32768 ! All integer bulk data (local indices, counts, maps); not ids
        integer, parameter :: EX_ALL_INT64_API =  8192+16384+32768 ! All of the above...

!
!  -- For EXSETOPT Calls
        integer, parameter :: EX_OPT_MAX_NAME_LENGTH = 1
        integer, parameter :: EX_OPT_COMPRESSION_TYPE = 2 !Currently not used. GZip by default
        integer, parameter :: EX_OPT_COMPRESSION_LEVEL = 3 ! 0 (disabled/fastest) ... 9 (best/slowest)
        integer, parameter :: EX_OPT_COMPRESSION_SHUFFLE = 4 ! 0 (disabled); 1 (enabled)
        integer, parameter :: EX_OPT_INTEGER_SIZE_API = 5 ! See *_INT64_* values above
        integer, parameter :: EX_OPT_INTEGER_SIZE_DB  = 6 ! (query only)
!
!  --INQUIRE FLAGS
        integer, parameter :: EX_INQ_FILE_TYPE                  = 1 ! EXODUS file type (deprecated) */
        integer, parameter :: EX_INQ_API_VERS                   = 2 ! API version number (float) */
        integer, parameter :: EX_INQ_DB_VERS                    = 3 ! database version number (float) */
        integer, parameter :: EX_INQ_TITLE                      = 4 ! database title. MAX_LINE_LENGTH+1 char* size */
        integer, parameter :: EX_INQ_DIM                        = 5 ! number of dimensions */
        integer, parameter :: EX_INQ_NODES                      = 6 ! number of nodes    */
        integer, parameter :: EX_INQ_ELEM                       = 7 ! number of elements */
        integer, parameter :: EX_INQ_ELEM_BLK                   = 8 ! number of element blocks */
        integer, parameter :: EX_INQ_NODE_SETS                  = 9 ! number of node sets*/
        integer, parameter :: EX_INQ_NS_NODE_LEN                = 10 ! length of node set node list */
        integer, parameter :: EX_INQ_SIDE_SETS                  = 11 ! number of side sets*/
        integer, parameter :: EX_INQ_SS_NODE_LEN                = 12 ! length of side set node list */
        integer, parameter :: EX_INQ_SS_ELEM_LEN                = 13 ! length of side set element list */
        integer, parameter :: EX_INQ_QA                         = 14 ! number of QA records */
        integer, parameter :: EX_INQ_INFO                       = 15 ! number of info records */
        integer, parameter :: EX_INQ_TIME                       = 16 ! number of time steps in the database */
        integer, parameter :: EX_INQ_EB_PROP                    = 17 ! number of element block properties */
        integer, parameter :: EX_INQ_NS_PROP                    = 18 ! number of node set properties */
        integer, parameter :: EX_INQ_SS_PROP                    = 19 ! number of side set properties */
        integer, parameter :: EX_INQ_NS_DF_LEN                  = 20 ! length of node set distribution factor list*/
        integer, parameter :: EX_INQ_SS_DF_LEN                  = 21 ! length of side set distribution factor list*/
        integer, parameter :: EX_INQ_LIB_VERS                   = 22 ! API Lib vers number (float) */
        integer, parameter :: EX_INQ_EM_PROP                    = 23 ! number of element map properties */
        integer, parameter :: EX_INQ_NM_PROP                    = 24 ! number of node map properties */
        integer, parameter :: EX_INQ_ELEM_MAP                   = 25 ! number of element maps */
        integer, parameter :: EX_INQ_NODE_MAP                   = 26 ! number of node maps*/
        integer, parameter :: EX_INQ_EDGE                       = 27 ! number of edges    */
        integer, parameter :: EX_INQ_EDGE_BLK                   = 28 ! number of edge blocks */
        integer, parameter :: EX_INQ_EDGE_SETS                  = 29 ! number of edge sets   */
        integer, parameter :: EX_INQ_ES_LEN                     = 30 ! length of concat edge set edge list       */
        integer, parameter :: EX_INQ_ES_DF_LEN                  = 31 ! length of concat edge set dist factor list*/
        integer, parameter :: EX_INQ_EDGE_PROP                  = 32 ! number of properties stored per edge block    */
        integer, parameter :: EX_INQ_ES_PROP                    = 33 ! number of properties stored per edge set      */
        integer, parameter :: EX_INQ_FACE                       = 34 ! number of faces */
        integer, parameter :: EX_INQ_FACE_BLK                   = 35 ! number of face blocks */
        integer, parameter :: EX_INQ_FACE_SETS                  = 36 ! number of face sets */
        integer, parameter :: EX_INQ_FS_LEN                     = 37 ! length of concat face set face list */
        integer, parameter :: EX_INQ_FS_DF_LEN                  = 38 ! length of concat face set dist factor list*/
        integer, parameter :: EX_INQ_FACE_PROP                  = 39 ! number of properties stored per face block */
        integer, parameter :: EX_INQ_FS_PROP                    = 40 ! number of properties stored per face set */
        integer, parameter :: EX_INQ_ELEM_SETS                  = 41 ! number of element sets */
        integer, parameter :: EX_INQ_ELS_LEN                    = 42 ! length of concat element set element list       */
        integer, parameter :: EX_INQ_ELS_DF_LEN                 = 43 ! length of concat element set dist factor list*/
        integer, parameter :: EX_INQ_ELS_PROP                   = 44 ! number of properties stored per elem set      */
        integer, parameter :: EX_INQ_EDGE_MAP                   = 45 ! number of edge maps                     */
        integer, parameter :: EX_INQ_FACE_MAP                   = 46 ! number of face maps                     */
        integer, parameter :: EX_INQ_COORD_FRAMES               = 47 ! number of coordinate frames */
        integer, parameter :: EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH = 48 ! size of MAX_NAME_LENGTH dimension on database */
        integer, parameter :: EX_INQ_DB_MAX_USED_NAME_LENGTH    = 49 ! size of MAX_NAME_LENGTH dimension on database */
        integer, parameter :: EX_INQ_MAX_READ_NAME_LENGTH       = 50 ! client-specified max size of returned names */
        integer, parameter :: EX_INQ_DB_FLOAT_SIZE    = 51 ! size of floating-point values stored on database */
        integer, parameter :: EX_INQ_NUM_CHILD_GROUPS = 52 ! number of groups contained in this (exoid) group */
        integer, parameter :: EX_INQ_GROUP_PARENT     = 53 ! id of parent of this (exoid) group; returns exoid if at root */
        integer, parameter :: EX_INQ_GROUP_ROOT = 54 ! id of root group "/" of this (exoid) group; returns exoid if at root */
        integer, parameter :: EX_INQ_GROUP_NAME_LEN = 55 ! length of name of group exoid */
        integer, parameter :: EX_INQ_GROUP_NAME = 56 ! name of group exoid. "/" returned for root group (char* GROUP_NAME_LEN+1 size) */
        integer, parameter :: EX_INQ_FULL_GROUP_NAME_LEN = 57 ! length of full path name of this (exoid) group */
        integer, parameter :: EX_INQ_FULL_GROUP_NAME     = 58 ! full "/"-separated path name of this (exoid) group */
        integer, parameter :: EX_INQ_THREADSAFE          = 59 ! Returns 1 if library is thread-safe; 0 otherwise */
        integer, parameter :: EX_INQ_ASSEMBLY            = 60 ! number of assemblies */
        integer, parameter :: EX_INQ_BLOB                = 61 ! number of blobs */
        integer, parameter :: EX_INQ_NUM_NODE_VAR        = 62 ! number of nodal variables */
        integer, parameter :: EX_INQ_NUM_EDGE_BLOCK_VAR  = 63 ! number of edge block variables */
        integer, parameter :: EX_INQ_NUM_FACE_BLOCK_VAR  = 64 ! number of face block variables */
        integer, parameter :: EX_INQ_NUM_ELEM_BLOCK_VAR  = 65 ! number of element block variables */
        integer, parameter :: EX_INQ_NUM_NODE_SET_VAR    = 66 ! number of node set variables */
        integer, parameter :: EX_INQ_NUM_EDGE_SET_VAR    = 67 ! number of edge set variables */
        integer, parameter :: EX_INQ_NUM_FACE_SET_VAR    = 68 ! number of face set variables */
        integer, parameter :: EX_INQ_NUM_ELEM_SET_VAR    = 69 ! number of element set variables */
        integer, parameter :: EX_INQ_NUM_SIDE_SET_VAR    = 70 ! number of sideset variables */
        integer, parameter :: EX_INQ_NUM_GLOBAL_VAR      = 71 ! number of global variables */
        integer, parameter :: EX_INQ_FILE_FORMAT         = 72 ! netCDF file format */
        integer, parameter :: EX_INQ_INVALID             = -1

        integer, parameter :: EXTYPE = 1  ! --  inquire EXODUS II file type
        integer, parameter :: EXVERS = 2  ! --  inquire API version number
        integer, parameter :: EXDBVR = 3  ! --  inquire database version number
        integer, parameter :: EXTITL = 4  ! --  inquire database title
        integer, parameter :: EXDIM  = 5  ! --  inquire number of dimensions
        integer, parameter :: EXNODE = 6  ! --  inquire number of nodes
        integer, parameter :: EXELEM = 7  ! --  inquire number of elements
        integer, parameter :: EXELBL = 8  ! --  inquire number of element blocks
        integer, parameter :: EXNODS = 9  ! --  inquire number of node sets
        integer, parameter :: EXNSNL = 10 ! --  inquire length of node set node list
        integer, parameter :: EXSIDS = 11 ! --  inquire number of side sets
        integer, parameter :: EXSSNL = 12 ! --  inquire length of side set node list
        integer, parameter :: EXSSEL = 13 ! --  inquire length of side set element list
        integer, parameter :: EXQA   = 14 ! --  inquire number of QA records
        integer, parameter :: EXINFO = 15 ! --  inquire number of info records
        integer, parameter :: EXTIMS = 16 ! --  inquire number of time steps in the database
        integer, parameter :: EXNEBP = 17 ! --  inquire number of element block properties
        integer, parameter :: EXNNSP = 18 ! --  inquire number of node set properties
        integer, parameter :: EXNSSP = 19 ! --  inquire number of side set properties
        integer, parameter :: EXNSDF = 20 ! --  inquire length of node set distribution factor list
        integer, parameter :: EXSSDF = 21 ! --  inquire length of side set distribution factor list
        integer, parameter :: EXLBVR = 22 ! --  inquire API Lib vers number
        integer, parameter :: EXNEMP = 23 ! --  inquire number of element map properties
        integer, parameter :: EXNNMP = 24 ! --  inquire number of node map properties
        integer, parameter :: EXNEM  = 25 ! --  inquire number of element maps
        integer, parameter :: EXNNM  = 26 ! --  inquire number of node maps
        integer, parameter :: EXNEDG = 27 ! --  inquire number of edges
        integer, parameter :: EXNEDB = 28 ! --  inquire number of edge blocks
        integer, parameter :: EXNEDS = 29 ! --  inquire number of edge sets
        integer, parameter :: EXESL  = 30 ! --  inquire length of concat edge set edge list
        integer, parameter :: EXESDL = 31 ! --  inquire length of concat edge set dist factor list
        integer, parameter :: EXEDBP = 32 ! --  inquire number of properties stored per edge block
        integer, parameter :: EXEDSP = 33 ! --  inquire number of properties stored per edge set
        integer, Parameter :: EXNFAC = 34 ! --  inquire number of faces
        integer, parameter :: EXNFB  = 35 ! --  inquire number of face blocks
        integer, parameter :: EXNFS  = 36 ! --  inquire number of face sets
        integer, parameter :: EXFSL  = 37 ! --  inquire length of concat face set face list
        integer, parameter :: EXFSDL = 38 ! --  inquire length of concat face set dist factor list
        integer, parameter :: EXFABP = 39 ! --  inquire number of properties stored per face block
        integer, parameter :: EXFASP = 40 ! --  inquire number of properties stored per face set
        integer, parameter :: EXNELS = 41 ! --  inquire number of element sets
        integer, parameter :: EXELSL = 42 ! --  inquire length of concat element set element list
        integer, parameter :: EXELSD = 43 ! --  inquire length of concat element set dist factor list
        integer, parameter :: EXELSP = 44 ! --  inquire number of properties stored per elem set
        integer, parameter :: EXNEDM = 45 ! --  inquire number of edge maps
        integer, parameter :: EXNFAM = 46 ! --  inquire number of face maps
        integer, parameter :: EXNCF  = 47 ! --  inquire number of coordinate frames
        integer, parameter :: EXDBMXALNM = 48 ! --  inquire size of MAX_NAME_LENGTH dimension on database
        integer, parameter :: EXDBMXUSNM = 49 ! --  inquire max size of any name existing on database
        integer, parameter :: EXMXRDNM = 50 ! --  inquire user-specified max size of names returned; truncated if longer
!
!         Define limits
!
!     MXSTLN - Maximum string length: used for store strings
        integer, parameter :: MXSTLN = 32
!     MXLNLN - Maximum line string length
        integer, parameter :: MXLNLN = 80
!     MXNAME - Maximum possible length of a name (due to netcdf NC_MAX_NAME)
        integer, parameter :: MXNAME = 256
!
!     properties
        integer, parameter :: EXEBLK = 1  ! EXEBLK - element block property code
        integer, parameter :: EXNSET = 2  ! EXNSPR - node set property code
        integer, parameter :: EXSSET = 3  ! EXSSPR - side set property code
        integer, parameter :: EXEMAP = 4  ! EXEMAP - element map property code
        integer, parameter :: EXNMAP = 5  ! EXNMAP - node map property code
!
!     Define values for type
        integer, parameter :: EX_NODAL       = 14 !     nodal "block" for variables
        integer, parameter :: EX_NODE_BLOCK  = 14 !     alias for EX_NODAL
        integer, parameter :: EX_NODE_SET    =  2 !     node set property code
        integer, parameter :: EX_EDGE_BLOCK  =  6 !     edge block property code
        integer, parameter :: EX_EDGE_SET    =  7 !     edge set property code
        integer, parameter :: EX_FACE_BLOCK  =  8 !     face block property code
        integer, parameter :: EX_FACE_SET    =  9 !     face set property code
        integer, parameter :: EX_ELEM_BLOCK  =  1 !     element block property code
        integer, parameter :: EX_ELEM_SET    = 10 !     face set property code
        integer, parameter :: EX_SIDE_SET    =  3 !     side set property code

        integer, parameter :: EX_ELEM_MAP    =  4 !     element map property code
        integer, parameter :: EX_NODE_MAP    =  5 !     node map property code
        integer, parameter :: EX_EDGE_MAP    = 11 !     edge map property code
        integer, parameter :: EX_FACE_MAP    = 12 !     face map property code

        integer, parameter :: EX_GLOBAL      = 13 !     global "block" for variables
        integer, parameter :: EX_COORDINATE  = 15 !     kluge so some internal wrapper functions work
        integer, parameter :: EX_INVALID     = -1
!
!     Define options codes
        integer, parameter :: EXVRBS = 1  ! EXVRBS - verbose mode message flag
        integer, parameter :: EXDEBG = 2  ! EXDEBG - debug mode message flag
        integer, parameter :: EXABRT = 4  ! EXABRT - abort mode flag
!
!     Define coordinate frames
        integer, parameter :: EXCFREC = 1 ! Rectangular
        integer, parameter :: EXCFCYL = 2 ! Cylindrical
        integer, parameter :: EXCFSPH = 3 ! Spherical
!
!     Define error return codes
        integer, parameter :: EXFATL = -1 ! EXFATL - fatal error flag def
        integer, parameter :: EXOK = 0    ! EXOK - no err (success) flag def
        integer, parameter :: EXWARN = 1  ! EXWARN - warning flag def
!
!     Define error return codes and flags
        integer, parameter :: EXMEMF =  1000 ! EXMEMF - memory allocation failure flag def
        integer, parameter :: EXBFMD =  1001 ! EXBFMD - wrong file mode
        integer, parameter :: EXBFID =  1002 ! EXBFID - bad file id def
        integer, parameter :: EXBTID =  1004 ! EXBTID - id table lookup failed
        integer, parameter :: EXBPRM =  1005 ! EXBPRM - bad parameter passed
        integer, parameter :: EXNULL = -1006 ! EXNULL - null entity found
        integer, parameter :: EXPMSG = -1000 ! EXPMSG - message print code - no error implied
        integer, parameter :: EXLMSG = -1001 ! EXLMSG - print last error message msg code

        integer, parameter :: lk = selected_int_kind(15)
        integer*8, parameter :: EX_FILE_ID_MASK = 4294901760_lk ! xFFFF0000
