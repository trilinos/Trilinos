C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C============================================================================
      SUBROUTINE STCLST(MODULE)
C============================================================================

C   --*** STCLST ***  (BLOT) set color list mapping

C --Parameters:

C --  MODULE - character name of the subroutine module setting the colors

C  PURPOSE:  SETS THE DESIRED MAPPING OF THE DEFAULT COLORS.  THE ROUTINE
C            WAS WRITTEN TO ALLOW THE USER TO CHANGE THE MAPPING FOR
C            DIFFERENT PROGRAMS, SPECIFICALLY MAPPING ONE WAY FOR CONTOURS
C            AND PAINTING AND ANOTHER FOR PLOTTING.

      CHARACTER*(*) MODULE
      include 'params.blk'
      include 'cmap-lst.blk'
      CHARACTER*1 LSTMOD
      SAVE LSTMOD
      LOGICAL PLTICL

      CHARACTER*(MXSTLN) LIST1(NCOLOR), LIST2(NCOLOR)

      DATA LSTMOD / ' '/
      DATA LIST1 / 'BLACK   ', 'WHITE   ', 'RED     ',
     &   'GREEN   ', 'YELLOW  ', 'BLUE    ',
     &   'CYAN    ', 'MAGENTA ', '        ' /
      DATA LIST2 / 'BLACK   ', 'WHITE   ', 'BLUE    ',
     &   'MAGENTA ', 'CYAN    ', 'GREEN   ',
     &   'YELLOW  ', 'RED     ', '        ' /

C SET COLLST TO EITHER LIST DEPENDING ON THE VALUE OF MODULE

      IF(LSTMOD .NE. MODULE(1:1)) THEN
         LSTMOD = MODULE(1:1)

         IF(MODULE(1:1) .EQ. 'D' .OR. MODULE(1:1) .EQ. 'd') THEN
            DO 10 I = 1, NCOLOR
               COLLST(I) = LIST2(I)
10          CONTINUE
         ELSE
            DO 20 I = 1, NCOLOR
               COLLST(I) = LIST1(I)
20          CONTINUE
         END IF

         DO 30 I = 1, NCOLOR-2
            IF (PLTICL (COLLST(I+2), RCOLOR)) THEN
               COLMAP(I) = RCOLOR
            ELSE
               COLMAP(I) = -1
            ENDIF
30       CONTINUE
      END IF

      RETURN
      END
