C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPSHOW (SHOTYP, NAMES, MAPEL, MAPND)
C=======================================================================

C   --*** TPSHOW *** (TPLOT) Display TPLOT parameter information
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --TPSHOW displays the TPLOT plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   TYPLOT   - the curves to be plotted for the plot set
C   --   XYPLOT   -
C   --   PLOT     -
C   --   HARDCOPY -
C   --   NEUTRAL  -
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMES - IN - the variable names
C   --
C   --Common Variables:
C   --   Uses NTPVAR, TIMPLT, ITVID, ITVNE of /TPVARS/

      include 'params.blk'
      include 'tpvars.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMES(*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL ISABRT
      CHARACTER*(1024) PV1, PV2
      CHARACTER*2 STR2

      IF ((SHOTYP .EQ. 'TYPLOT') .OR. (SHOTYP .EQ. 'XYPLOT')
     &   .OR. (SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')
     &   .OR. (SHOTYP .EQ. 'NEUTRAL')) THEN
         N = 1
         DO 100 NP = 1, NTPCRV
            IF (ISABRT ()) RETURN
            IF (TIMPLT) THEN
               PV1 = 'TIME'
            ELSE
               CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N),
     &            PV1, MAPEL, MAPND)
               N = N + 1
            END IF
            CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *        MAPEL, MAPND)
            N = N + 1
            WRITE (STR2, '(I2)', IOSTAT=IDUM) NP
            WRITE (*, 10000) 'Curve ', STR2, ' :  ', PV2(:LENSTR(PV2)),
     &         '  -vs-  ', PV1(:LENSTR(PV1))
  100    CONTINUE

      END IF

      RETURN

10000  FORMAT (1X, 10A)
      END
