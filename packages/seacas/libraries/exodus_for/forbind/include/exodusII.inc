! -*- mode: fortran -*-
!
! Copyright 1994 Sandia Corporation. Under the terms of Contract
! DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
! retains certain rights in this software.
!
! See the COPYRIGHT file for copying and redistribution conditions.
!
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
        integer excre,exopen,exlgmd
        external excre, exopen, exlgmd
!
!       Define access modes
!
        integer EXNOCL,EXCLOB,EXREAD,EXWRIT,EXWRT,EXNORM,EXLARG
        integer EXNET4, EXNOSH, EXSHAR

        parameter (EXNOCL = 0)   !       EXNOCL - File create with no overwrite (no clobber), default
        parameter (EXCLOB = 1)   !       EXCLOB - File create with overwrite (clobber)
        parameter (EXNORM = 2)   !       EXNORM - Disable large model modifications
        parameter (EXLARG = 4)   !       EXLARG - Large model modifications
        parameter (EXNET4 = 8)   !       EXNET4 - use the hdf5-based netcdf4 output
        parameter (EXNOSH = 16)  !       EXNOSH - Do not open netcdf file in "share" mode
        parameter (EXSHAR = 32)  !       EXSHAR - Do open netcdf file in "share" mode

        parameter (EXREAD = 0)   !       EXREAD - File open for read
        parameter (EXWRIT = 1)   !       EXWRIT - File open for write
        parameter (EXWRT  = 1)   !       EXWRT  - File open for write (alternate to EXWRIT)
!
!  --INQUIRE FLAGS
      integer EXTYPE, EXVERS, EXDBVR, EXTITL, EXDIM, EXNODE,            &
     &        EXELEM, EXELBL, EXNODS, EXSIDS, EXQA,  EXINFO,            &
     &        EXTIMS, EXSSNL, EXSSEL, EXNSNL, EXNEBP, EXNNSP,           &
     &        EXNSSP, EXNSDF, EXSSDF, EXLBVR, EXNEMP, EXNNMP,           &
     &        EXNEM,  EXNNM,  EXNEDG, EXNEDB, EXNEDS, EXESL ,           &
     &        EXESDL, EXEDBP, EXEDSP, EXNFAC, EXNFB , EXNFS ,           &
     &        EXFSL , EXFSDL, EXFABP, EXFASP, EXNELS, EXELSL,           &
     &        EXELSD, EXELSP, EXNEDM, EXNFAM, EXNCF,  EXDBMXALNM,       &
     &        EXDBMXUSNM, EXMXRDNM

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
        integer MXSTLN, MXLNLN
!       MXSTLN - Maximum string length: used for store strings
        parameter (MXSTLN = 32)
!       MXLNLN - Maximum line string length
        parameter (MXLNLN = 80)
!
!       properties
!
        integer EXEBLK, EXNSET, EXSSET, EXEMAP, EXNMAP

        parameter (EXEBLK = 1)   ! EXEBLK - element block property code
        parameter (EXNSET = 2)   ! EXNSPR - node set property code
        parameter (EXSSET = 3)   ! EXSSPR - side set property code
        parameter (EXEMAP = 4)   ! EXEMAP - element map property code
        parameter (EXNMAP = 5)   ! EXNMAP - node map property code
!
!       Define options codes
!
        integer EXVRBS, EXDEBG, EXABRT

        parameter (EXVRBS = 1)   ! EXVRBS - verbose mode message flag
        parameter (EXDEBG = 2)   ! EXDEBG - debug mode message flag
        parameter (EXABRT = 4)   ! EXABRT - abort mode flag
!
!       Define coordinate frames
!
        integer EXCFREC, EXCFCYL, EXCFSPH
        
        parameter (EXCFREC = 1)     ! Rectangular
        parameter (EXCFCYL = 2)     ! Cylindrical
        parameter (EXCFSPH = 3)     ! Spherical

!
!       Define error return codes
!
        integer EXFATL, EXOK, EXWARN

        parameter (EXFATL = -1)  ! EXFATL - fatal error flag def
        parameter (EXOK = 0)     ! EXOK - no err (success) flag def
        parameter (EXWARN = 1)   ! EXWARN - warning flag def
!
!       Define error return codes and flags
!
        integer EXMEMF, EXBFMD, EXBFID, EXBTID, EXBPRM, EXNULL,         &
     &          EXPMSG, EXLMSG

        parameter (EXMEMF =  1000)   ! EXMEMF - memory allocation failure flag def
        parameter (EXBFMD =  1001)   ! EXBFMD - wrong file mode
        parameter (EXBFID =  1002)   ! EXBFID - bad file id def
        parameter (EXBTID =  1004)   ! EXBTID - id table lookup failed
        parameter (EXBPRM =  1005)   ! EXBPRM - bad parameter passed
        parameter (EXNULL = -1006)   ! EXNULL - null entity found
        parameter (EXPMSG = -1000)   ! EXPMSG - message print code - no error implied
        parameter (EXLMSG = -1001)   ! EXLMSG - print last error message msg code
