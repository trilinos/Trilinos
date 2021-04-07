C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDCORD (NDB, NDIM, NUMNP, CORD, NAMECO, ISEOF, NAMLEN)
C=======================================================================

C   --*** RDCORD *** (EXPLORE) Read database coordinates
C   --
C   --RDCORD reads the coordinate array from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   CORD - OUT - the coordinates
C   --   ISEOF - IN/OUT - set true if end of file read
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      include 'exodusII.inc'

      REAL CORD(NUMNP,*)
      CHARACTER*(NAMLEN) NAMECO(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIREA (NDIM*NUMNP, 0.0, CORD)
      CALL INISTR (NDIM, '--------------------------------', NAMECO)

      if (NDIM .gt. 0) then
         call exgcor(ndb, cord(1,1), cord(1,2), cord(1,3), ierr)
         if (ierr .ne. 0) go to 100

         call exgcon(ndb, nameco, ierr)
         if (ierr .ne. 0) go to 105
      end if

      RETURN

 100  CONTINUE
      WRITE (ERRMSG, 10000) 'COORDINATES'
      GOTO 110
 105  CONTINUE
      WRITE (ERRMSG, 10000) 'COORDINATE NAMES'
      GOTO 110
 110  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
      RETURN

10000 FORMAT (A)
      END

