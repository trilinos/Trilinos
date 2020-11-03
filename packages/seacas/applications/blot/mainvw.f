C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION MAINVW ()
C=======================================================================

C   --*** MAINVW *** (MESH) Return the "main" view number
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --MAINVW returns the number of the "main" view for symmetric views.
C   --For non-symmetric views, zero is returned.
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/
C   --   Uses XISSYM, YISSYM, LFTSYM, BOTSYM of /VIEWS/

      include 'mshopt.blk'
      include 'views.blk'

      IF (XISSYM .OR. YISSYM) THEN
         IF ((.NOT. XISSYM) .OR. LFTSYM) THEN
            IF (BOTSYM) THEN
               MAINVW = 2
            ELSE
               MAINVW = 4
            END IF
         ELSE IF ((.NOT. YISSYM) .OR. BOTSYM) THEN
            IF (LFTSYM) THEN
               MAINVW = 2
            ELSE
               MAINVW = 1
            END IF
         ELSE
            MAINVW = 3
         END IF
      ELSE
         MAINVW = 0
      END IF

      RETURN
      END
