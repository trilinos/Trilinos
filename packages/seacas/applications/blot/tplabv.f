C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPLABV (NPLT, IVAR, NAME, NE, LABSTR, MAPEL, MAPND)
C=======================================================================

C   --*** TPLABV *** (TPLOT) Get a plot label
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --TPLABV returns a plot label.
C   --
C   --Parameters:
C   --   NPLT - IN - the label type:
C   --      -1 = name and number (long form)
C   --       0 = name and number (short form)
C   --   IVAR - IN - the variable index which determines type
C   --   NAME - IN - the variable name
C   --   NE - IN - the variable number
C   --   LABSTR - OUT - the label string
C   --
C   --Common Variables:
C   --   Uses NVARNP, NVAREL of /DBNUMS/

      CHARACTER*(*) NAME
      CHARACTER*(*) LABSTR
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER TYP

      CALL DBVTYP_BL (IVAR, TYP, IDUM)

      IF (NPLT .LE. -1) THEN
         IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME
         ELSE IF (TYP .EQ. 'N') THEN
           WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'at NODE',
     *       MAPND(NE)
         ELSE IF (TYP .EQ. 'E') THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'at ELEMENT',
     *       MAPEL(NE)
         END IF

      ELSE

         IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME
         ELSE IF (TYP .EQ. 'N') THEN
            WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'NODE',
     *       MAPND(NE)
         ELSE IF (TYP .EQ. 'E') THEN
           WRITE (LABSTR, 10000, IOSTAT=IDUM) NAME, 'ELEM',
     *       MAPEL(NE)
         END IF
      END IF

      CALL SQZSTR (LABSTR, LSTR)

      RETURN
10000  FORMAT (A, :, ' ', A, I12)
      END
