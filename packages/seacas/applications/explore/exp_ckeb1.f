C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CKEB1 (IEL0, IELB, IDELB, NUMELB, NUMLNK, NUMNP, LINK,
     *  NODUSE)
C=======================================================================

C   --*** CKEB1 *** (EXPLORE) Check database element block connectivity
C   --
C   --CKEB1 checks that the database element block connectivity is within
C   --the nodal range.
C   --
C   --Parameters:
C   --   IEL0 - IN - the number of the first element in this block - 1
C   --   IELB - IN - the number of this element block
C   --   IDELB - IN - the element block ID for this block
C   --   NUMELB - IN - the number of elements for this block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMNP - IN - the number of nodes
C   --   LINK - IN - the connectivity array for this block
C   --   NODUSE - OUT - scratch array to determine whether all nodes used by an element.

      include 'exp_errcnt.blk'
      INTEGER LINK(NUMLNK,*)
      INTEGER NODUSE(*)
      CHARACTER*132 STRA

      IF (NUMELB .LE. 0) GOTO 110
      IF (NUMLNK .LE. 0) GOTO 110

      NERR = 0

      DO 100 NE = 1, NUMELB
         CALL CHKRNG (LINK(1,NE), NUMLNK, NUMNP, NZERO, IERR)
         IF (IERR .GT. 0 .OR. NZERO .GT. 0) THEN
            IF (NERR .EQ. 0) THEN
               WRITE (*, 10000, IOSTAT=IDUM)
     &            'Connectivity Problems', IELB, IDELB
            END IF
            if (nerr .lt. maxerrs .or. maxerrs .le. 0) then
               WRITE (*, 10010, IOSTAT=IDUM)
     &              NE, NE+IEL0, (LINK(I,NE), I=1,NUMLNK)
            else if (nerr .eq. maxerrs .and. maxerrs .gt. 0) then
               call prterr('CMDSPEC',
     $              '...skipping additional errors...')
            end if
            NERR = NERR + 1
         END IF

         DO 90 I=1,NUMLNK
           NODE = LINK(I,NE)
           IF (NODE .LE. NUMNP) THEN
C ... Note that if NODE is out of range, the error message above will have
C     already been printed, so we don't print anything here.
             NODUSE(NODE) = 1
           END IF
 90      CONTINUE

 100  CONTINUE
      if (nerr .gt. 0) then
         write (stra, 10020) nerr, idelb
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
      end if

  110 CONTINUE
      RETURN

10000  FORMAT (/, 1X, '     #  elem      ', A, ' for block #', I12,
     &   ',', I12, ' = ID')
10010  FORMAT (1X, I12, I12, 5X, 6I12, :, /,
     &   (26X, 6I12))
10020    FORMAT('ELEMENT CONNECTIVITY ERROR: Found ',I12,
     $      ' errors in element connectivity check for element block '
     $      , i10)
      END
