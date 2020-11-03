C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE SELSSETN (NUMSEL, IXSEL, LISET, IDSSET,
     &      VALNAM, ISCRN, IA)
C=======================================================================

      include 'exodusII.inc'
      include 'exp_dbase.blk'
      include 'exp_dbnums.blk'

      INTEGER IXSEL(*)
      INTEGER LISET(0:*)
      INTEGER IDSSET(*)
      INTEGER ISCRN(*)
      INTEGER IA(*)

      CHARACTER*(*) VALNAM
      CHARACTER*40 STRA
      CHARACTER*132 MSG

      do i=1, numnp
        iscrn(i) = 0
      end do

      NUMSEL = 0
      DO IX = 1, LISET(0)
         INSS = LISET(IX)
         IDSS = IDSSET(INSS)

C     ... Get count of nodes in the sideset nodelist (includes duplicates)
         CALL EXGSNL(NDB, IDSS, NNESS, IERR)

         CALL MDRSRV('ISSCR', KNODSS, NNESS)
C ... This is too large, but easier to just guess than get value now.
         CALL MDRSRV('ISSCE', KNODSE, NNESS/3)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 280

         CALL EXGSSN(NDB, IDSS, IA(KNODSE), IA(KNODSS), IERR)

         DO I=0, NNESS-1
            NODE = IA(KNODSS + I)
            ISCRN(NODE) = 1
         end do

         call mddel('ISSCE')
         call mddel('ISSCR')
      end do

C ... Iterate 'iscrn' and select any node with iscrn(i) = 1
      do i = 1, numnp
         if (iscrn(i) .eq. 1) then
            numsel = numsel + 1
            ixsel(numsel) = i
         end if
      end do

      write (stra, 10000) numsel
      call pckstr(1, stra)
      MSG = STRA(:lenstr(stra)) // ' ' // VALNAM // ' selected'
      call prterr('CMDSPEC', MSG)
10000 format(I12)
 280  continue
      RETURN
      END
