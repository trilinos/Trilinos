!
!  
!  (C) 1993 by Argonne National Laboratory and Mississipi State University.
!      All rights reserved.  See COPYRIGHT in top-level directory.
!
!
! user include file for MPI programs, with no dependencies 
!
! It really is not possible to make a perfect include file that can
! be used by both F77 and F90 compilers, but this is close.  We have removed
! continuation lines (allows free form input in F90); systems whose
! Fortran compilers support ! instead of just C or * for comments can
! globally replace a C in the first column with !; the resulting file
! should work for both Fortran 77 and Fortran 90.
!
! If your Fortran compiler supports ! for comments, you can run this 
! through sed with
!     sed -e 's/^C/\!/g'
!
! We have also removed the use of contractions (involving the single quote)
! character because some users use .F instead of .f files (to invoke the
! cpp preprocessor) and further, their preprocessor is determined to find
! matching single quote pairs (and probably double quotes; given the
! different rules in C and Fortran, this sounds like a disaster).  Rather than
! take the position that the poor users should get a better system, we
! have removed the text that caused problems.  Of course, the users SHOULD
! get a better system...
!
! return codes 
      INTEGER MPI_SUCCESS,MPI_ERR_BUFFER,MPI_ERR_COUNT,MPI_ERR_TYPE
      INTEGER MPI_ERR_TAG,MPI_ERR_COMM,MPI_ERR_RANK,MPI_ERR_ROOT
      INTEGER MPI_ERR_GROUP
      INTEGER MPI_ERR_OP,MPI_ERR_TOPOLOGY,MPI_ERR_DIMS,MPI_ERR_ARG
      INTEGER MPI_ERR_UNKNOWN,MPI_ERR_TRUNCATE,MPI_ERR_OTHER
      INTEGER MPI_ERR_INTERN,MPI_ERR_IN_STATUS,MPI_ERR_PENDING
      INTEGER MPI_ERR_REQUEST, MPI_ERR_LASTCODE
      PARAMETER (MPI_SUCCESS=0,MPI_ERR_BUFFER=1,MPI_ERR_COUNT=2)
      PARAMETER (MPI_ERR_TYPE=3,MPI_ERR_TAG=4,MPI_ERR_COMM=5)
      PARAMETER (MPI_ERR_RANK=6,MPI_ERR_ROOT=7,MPI_ERR_GROUP=8)
      PARAMETER (MPI_ERR_OP=9,MPI_ERR_TOPOLOGY=10,MPI_ERR_DIMS=11)
      PARAMETER (MPI_ERR_ARG=12,MPI_ERR_UNKNOWN=13)
      PARAMETER (MPI_ERR_TRUNCATE=14,MPI_ERR_OTHER=15)
      PARAMETER (MPI_ERR_INTERN=16,MPI_ERR_IN_STATUS=17)
      PARAMETER (MPI_ERR_PENDING=18,MPI_ERR_REQUEST=19)
      PARAMETER (MPI_ERR_LASTCODE=4114)
!
      INTEGER MPI_UNDEFINED
      parameter (MPI_UNDEFINED = (-32766))
!
      INTEGER MPI_GRAPH, MPI_CART
      PARAMETER (MPI_GRAPH = 1, MPI_CART = 2)
      INTEGER  MPI_PROC_NULL
      PARAMETER ( MPI_PROC_NULL = (-1) )
!
      INTEGER MPI_BSEND_OVERHEAD
      PARAMETER ( MPI_BSEND_OVERHEAD = 512 )

      INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
      PARAMETER(MPI_SOURCE=2, MPI_TAG=3, MPI_ERROR=4)
      INTEGER MPI_STATUS_SIZE
      PARAMETER (MPI_STATUS_SIZE=5)
      INTEGER MPI_MAX_PROCESSOR_NAME, MPI_MAX_ERROR_STRING
      PARAMETER (MPI_MAX_PROCESSOR_NAME=256)
      PARAMETER (MPI_MAX_ERROR_STRING=512)
      INTEGER MPI_MAX_NAME_STRING
      PARAMETER (MPI_MAX_NAME_STRING=63)
!
      INTEGER MPI_COMM_NULL
      PARAMETER (MPI_COMM_NULL=0)
!
      INTEGER MPI_DATATYPE_NULL
      PARAMETER (MPI_DATATYPE_NULL = 0)
      
      INTEGER MPI_ERRHANDLER_NULL
      PARAMETER (MPI_ERRHANDLER_NULL = 0)
      
      INTEGER MPI_GROUP_NULL
      PARAMETER (MPI_GROUP_NULL = 0)
      
      INTEGER MPI_KEYVAL_INVALID
      PARAMETER (MPI_KEYVAL_INVALID = 0)
      
      INTEGER MPI_REQUEST_NULL
      PARAMETER (MPI_REQUEST_NULL = 0)
! 
      INTEGER MPI_IDENT, MPI_CONGRUENT, MPI_SIMILAR, MPI_UNEQUAL
      PARAMETER (MPI_IDENT=0, MPI_CONGRUENT=1, MPI_SIMILAR=2)
      PARAMETER (MPI_UNEQUAL=3)
!
!     MPI_BOTTOM needs to be a known address; here we put it at the
!     beginning of the common block.  The point-to-point and collective
!     routines know about MPI_BOTTOM, but MPI_TYPE_STRUCT as yet does not.
!
!     MPI_STATUS_IGNORE and MPI_STATUSES_IGNORE are similar objects
!
!     The types MPI_INTEGER1,2,4 and MPI_REAL4,8 are OPTIONAL.
!     Their values are zero if they are not available.  Note that
!     using these reduces the portability of code (though may enhance
!     portability between Crays and other systems)
!
      INTEGER MPI_TAG_UB, MPI_HOST, MPI_IO
      INTEGER MPI_BOTTOM, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
      INTEGER MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION 
      INTEGER MPI_COMPLEX, MPI_DOUBLE_COMPLEX,MPI_LOGICAL
      INTEGER MPI_CHARACTER, MPI_BYTE, MPI_2INTEGER, MPI_2REAL
      INTEGER MPI_2DOUBLE_PRECISION, MPI_2COMPLEX, MPI_2DOUBLE_COMPLEX
      INTEGER MPI_UB, MPI_LB
      INTEGER MPI_PACKED, MPI_WTIME_IS_GLOBAL
      INTEGER MPI_COMM_WORLD, MPI_COMM_SELF, MPI_GROUP_EMPTY
      INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD, MPI_LAND, MPI_BAND
      INTEGER MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC
      INTEGER MPI_MAXLOC
      INTEGER MPI_OP_NULL
      INTEGER MPI_ERRORS_ARE_FATAL, MPI_ERRORS_RETURN
!
      PARAMETER (MPI_ERRORS_ARE_FATAL=119)
      PARAMETER (MPI_ERRORS_RETURN=120)
!
      PARAMETER (MPI_COMPLEX=23,MPI_DOUBLE_COMPLEX=24,MPI_LOGICAL=25)
      PARAMETER (MPI_REAL=26,MPI_DOUBLE_PRECISION=27,MPI_INTEGER=28)
      PARAMETER (MPI_2INTEGER=29,MPI_2COMPLEX=30,MPI_2DOUBLE_COMPLEX=31)
      PARAMETER (MPI_2REAL=32,MPI_2DOUBLE_PRECISION=33,MPI_CHARACTER=1)
      PARAMETER (MPI_BYTE=3,MPI_UB=16,MPI_LB=15,MPI_PACKED=14)

      INTEGER MPI_ORDER_C, MPI_ORDER_FORTRAN 
      PARAMETER (MPI_ORDER_C=56, MPI_ORDER_FORTRAN=57)
      INTEGER MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_CYCLIC
      INTEGER MPI_DISTRIBUTE_NONE, MPI_DISTRIBUTE_DFLT_DARG
      PARAMETER (MPI_DISTRIBUTE_BLOCK=121, MPI_DISTRIBUTE_CYCLIC=122)
      PARAMETER (MPI_DISTRIBUTE_NONE=123)
      PARAMETER (MPI_DISTRIBUTE_DFLT_DARG=-49767)
      INTEGER MPI_MAX_INFO_KEY, MPI_MAX_INFO_VAL
      PARAMETER (MPI_MAX_INFO_KEY=255, MPI_MAX_INFO_VAL=1024)
      INTEGER MPI_INFO_NULL
      PARAMETER (MPI_INFO_NULL=0)

!
! Optional Fortran Types.  Configure attempts to determine these.  
!
      INTEGER MPI_INTEGER1, MPI_INTEGER2, MPI_INTEGER4, MPI_INTEGER8
      INTEGER MPI_INTEGER16
      INTEGER MPI_REAL4, MPI_REAL8, MPI_REAL16
      INTEGER MPI_COMPLEX8, MPI_COMPLEX16, MPI_COMPLEX32
      PARAMETER (MPI_INTEGER1=1,MPI_INTEGER2=4)
      PARAMETER (MPI_INTEGER4=6)
      PARAMETER (MPI_INTEGER8=13)
      PARAMETER (MPI_INTEGER16=0)
      PARAMETER (MPI_REAL4=10)
      PARAMETER (MPI_REAL8=11)
      PARAMETER (MPI_REAL16=0)
      PARAMETER (MPI_COMPLEX8=23)
      PARAMETER (MPI_COMPLEX16=24)
      PARAMETER (MPI_COMPLEX32=0)

      COMMON /MPIPRIV/ MPI_BOTTOM,MPI_STATUS_IGNORE,MPI_STATUSES_IGNORE      
!
!     Without this save, some Fortran implementations may make the common
!     dynamic!
!    
!     For a Fortran90 module, we might replace /MPIPRIV/ with a simple
!     SAVE MPI_BOTTOM
!
      SAVE /MPIPRIV/

      PARAMETER (MPI_MAX=100,MPI_MIN=101,MPI_SUM=102,MPI_PROD=103)
      PARAMETER (MPI_LAND=104,MPI_BAND=105,MPI_LOR=106,MPI_BOR=107)
      PARAMETER (MPI_LXOR=108,MPI_BXOR=109,MPI_MINLOC=110)
      PARAMETER (MPI_MAXLOC=111, MPI_OP_NULL=0)
!
      PARAMETER (MPI_GROUP_EMPTY=90,MPI_COMM_WORLD=91,MPI_COMM_SELF=92)
      PARAMETER (MPI_TAG_UB=80,MPI_HOST=82,MPI_IO=84)
      PARAMETER (MPI_WTIME_IS_GLOBAL=86)
!
      INTEGER MPI_ANY_SOURCE
      PARAMETER (MPI_ANY_SOURCE = (-2))
      INTEGER MPI_ANY_TAG
      PARAMETER (MPI_ANY_TAG = (-1))
!
      INTEGER MPI_VERSION, MPI_SUBVERSION
      PARAMETER (MPI_VERSION    = 1, MPI_SUBVERSION = 1)
!
!     All other MPI routines are subroutines
!     This may cause some Fortran compilers to complain about defined and
!     not used.  Such compilers should be improved.
!
      DOUBLE PRECISION MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
      EXTERNAL MPI_WTIME, MPI_WTICK, PMPI_WTIME, PMPI_WTICK
!
!     The attribute copy/delete subroutines are symbols that can be passed
!     to MPI routines
!
      EXTERNAL MPI_NULL_COPY_FN, MPI_NULL_DELETE_FN, MPI_DUP_FN
! 
!     $Id$    
! 
!     Copyright (C) 1997 University of Chicago. 
!     See COPYRIGHT notice in top-level directory.
!
! 
!    user include file for Fortran MPI-IO programs 
!
      INTEGER MPI_MODE_RDONLY, MPI_MODE_RDWR, MPI_MODE_WRONLY
      INTEGER MPI_MODE_DELETE_ON_CLOSE, MPI_MODE_UNIQUE_OPEN
      INTEGER MPI_MODE_CREATE, MPI_MODE_EXCL
      INTEGER MPI_MODE_APPEND
      PARAMETER (MPI_MODE_RDONLY=2, MPI_MODE_RDWR=8, MPI_MODE_WRONLY=4)
      PARAMETER (MPI_MODE_CREATE=1, MPI_MODE_DELETE_ON_CLOSE=16)
      PARAMETER (MPI_MODE_UNIQUE_OPEN=32, MPI_MODE_EXCL=64)
      PARAMETER (MPI_MODE_APPEND=128)
!
      INTEGER MPI_FILE_NULL
      PARAMETER (MPI_FILE_NULL=0)
!
      INTEGER MPI_MAX_DATAREP_STRING
      PARAMETER (MPI_MAX_DATAREP_STRING=128)
!
      INTEGER MPI_SEEK_SET, MPI_SEEK_CUR, MPI_SEEK_END
      PARAMETER (MPI_SEEK_SET=600, MPI_SEEK_CUR=602, MPI_SEEK_END=604)
!
      INTEGER MPIO_REQUEST_NULL
      PARAMETER (MPIO_REQUEST_NULL=0)
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
