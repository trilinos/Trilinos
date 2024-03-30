C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SETMOD (IVIEW, MMOD, MTYP)
C=======================================================================

C   --*** SETMOD *** (DETOUR) Set display mode and type
C   --   Written by Amy Gilkey - revised 04/06/88
C   --
C   --SETMOD sets the display mode and mode type of one or all views.
C   --Counts the number of variables needed, if changed
C   --
C   --Parameters:
C   --   IVIEW - IN - the view to be set, 0 for all
C   --   MMOD - IN - the display mode to be set
C   --   MTYP - IN - the display mode type to be set
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/
C   --   Sets MODDET, MODTYP, NNDVAR, NEDVAR of /DETOPT/

      include 'mshopt.blk'
      include 'detopt.blk'

      CHARACTER*(*) MMOD, MTYP

      IF (IVIEW .GE. 1) THEN
         IF (MMOD .EQ. 'NONE') THEN
            MODDET(IVIEW) = 'NONE'
            MODTYP(IVIEW) = ' '
         ELSE
            MODDET(IVIEW) = MMOD
            MODTYP(IVIEW) = MTYP
         END IF
      ELSE
         IF (MMOD .EQ. 'NONE') THEN
            DO 100 I = 1, 4
               MODDET(I) = 'NONE'
               MODTYP(I) = ' '
  100       CONTINUE
         ELSE
            DO 110 I = 1, 4
               IF ((MSHDEF(I) .EQ. 'NONE')
     &            .OR. (MSHDEF(I) .EQ. 'EMPTY')) THEN
                  MODDET(I) = 'NONE'
                  MODTYP(I) = ' '
               ELSE
                  MODDET(I) = MMOD
                  MODTYP(I) = MTYP
               END IF
  110       CONTINUE
         END IF
      END IF

C   --Calculate the number of variables needed
      CALL CNTVAR (MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR)

      RETURN
      END
