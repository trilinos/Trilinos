C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNTVAR (MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR)
C=======================================================================

C   --*** CNTVAR *** (DETOUR) Count the number of variables needed
C   --   Written by Amy Gilkey - revised 03/03/88
C   --
C   --CNTVAR counts the number of nodal and element database variables
C   --needed for the display modes.
C   --
C   --Parameters:
C   --   MODDET - IN - the modes for all views (as in /DETOPT/)
C   --   MODTYP - IN - the mode types for all views (as in /DETOPT/)
C   --   IDTVAR - IN - the current variables
C   --   NNDVAR - OUT - the number of nodal variables needed
C   --   NEDVAR - OUT - the number of element variables needed
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) MODDET(4), MODTYP(4)
      INTEGER IDTVAR(4)

      INTEGER NDEFVW, IXVW
      CHARACTER TYP

      NNDVAR = 0
      NEDVAR = 0
      DO 100 IVW = 1, NDEFVW (.FALSE.)
         IVIEW = IXVW (.FALSE., IVW)
         IF (MODDET(IVIEW) .EQ. 'CONTOUR') THEN
            NNDVAR = MAX (NNDVAR, 1)
            CALL DBVTYP_BL (IDTVAR(1), TYP, IDUM)
            IF (TYP .EQ. 'E') NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'ELEMCONT') THEN
            NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'VECTOR') THEN
            IF (MODTYP(IVIEW) .EQ. 'NODE') THEN
               NNDVAR = MAX (NNDVAR, NDIM)
            ELSE IF (MODTYP(IVIEW) .EQ. 'ELEMENT') THEN
               NEDVAR = MAX (NEDVAR, NDIM)
            ELSE IF ((MODTYP(IVIEW) .EQ. 'SIGMAX')
     &         .OR. (MODTYP(IVIEW) .EQ. 'SIGMIN')) THEN
               NEDVAR = MAX (NEDVAR, 3)
            END IF
         ELSE IF (MODDET(IVIEW) .EQ. 'SYMBOL') THEN
            NEDVAR = MAX (NEDVAR, 1)
         ELSE IF (MODDET(IVIEW) .EQ. 'GAUSS') THEN
            NEDVAR = MAX (NEDVAR, 4)
         END IF
  100 CONTINUE

      RETURN
      END
