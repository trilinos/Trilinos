C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNSHOW (SHOTYP, NAMES)
C=======================================================================

C   --*** LNSHOW *** (PATHLN) Display PATHLINE parameter information
C   --   Written by Amy Gilkey - revised 05/20/88
C   --
C   --LNSHOW displays the PATHLINE plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   LOCATION - the pathlines to be plotted for the plot set
C   --   PLOT     -
C   --   HARDCOPY -
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMES - IN - the variable names
C   --
C   --Common Variables:
C   --   Uses NLNCRV, ILVNE, ILVID of /LNVARS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'lnvars.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMES(*)

      LOGICAL ISABRT
      CHARACTER*(MXNAME) NAM(3)
      CHARACTER TYP
      CHARACTER*2 STR2
      CHARACTER*8 STRA
      CHARACTER*80 STRING

      IF ((SHOTYP .EQ. 'LOCATION')
     &   .OR. (SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')) THEN
         DO 110 NP = 1, NLNCRV
            IF (ISABRT ()) RETURN
            DO 100 IXY = 1, NDIM
               NAM(IXY) = NAMES(ILVID(IXY,NP))
  100       CONTINUE
            CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
            WRITE (STR2, '(I2)', IOSTAT=IDUM) NP
            WRITE (STRING, '(10A)') 'Pathline ', STR2, ' : ',
     &         (' ', NAM(I), I=1,NDIM), '^'
            LSTR = LENSTR (STRING) - 1
            IF ((TYP .EQ. 'N') .OR. (TYP .EQ. 'E')) THEN
               CALL INTSTR (1, 0, ILVNE(NP), STRA, L)
               IF (TYP .EQ. 'N') THEN
                  STRING(LSTR+1:) = ' Node ' // STRA(:L)
               ELSE
                  STRING(LSTR+1:) = ' Element ' // STRA(:L)
               END IF
            ELSE
               STRING(LSTR+1:LSTR+1) = ' '
            END IF
            WRITE (*, 10000) STRING(:LENSTR(STRING))
  110    CONTINUE

      END IF

      RETURN

10000  FORMAT (1X, 10A)
      END
