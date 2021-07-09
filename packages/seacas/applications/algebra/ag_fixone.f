C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE FIXONE (MAXNE, VARVAL)
C=======================================================================

C   --*** FIXONE *** (ALGEBRA) Handle first time step variables
C   --   Written by Amy Gilkey - revised 12/10/87
C   --
C   --FIXONE is called after the first time step is read in.  It copies
C   --the current variables to the last and first time step variables
C   --(if needed).  It also corrects the current variable pointer if
C   --only the first time step for a variable is needed.
C   --
C   --Parameters:
C   --   MAXNE  - IN     - the VARVAL dimension
C   --   VARVAL - IN/OUT - the input data and copied data
C   --
C   --Common Variables:
C   --   Uses ISTVAR of /VAR../
C   --   Uses NVARNP, NVAREL of /DBNUMS/

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      include 'ag_namlen.blk'
      include 'ag_var.blk'

      REAL VARVAL(MAXNE,*)

      DO 100 IVAR = 1, NUMINP
            IF (ISTVAR(IONETM,IVAR) .GT. 0) THEN
               IF (ISTVAR(ICURTM,IVAR) .LT. 0) THEN
                  ISTVAR(ICURTM,IVAR) = 0
               ELSE
                  IFROM = ISTVAR(ICURTM,IVAR)
                  ITO   = ISTVAR(IONETM,IVAR)
                  CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                         VARVAL(1,IFROM), VARVAL(1,ITO))
               END IF
            END IF
            IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
               IFROM = ISTVAR(ICURTM,IVAR)
               ITO = ISTVAR(ILSTTM,IVAR)
               CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                      VARVAL(1,IFROM), VARVAL(1,ITO))
            END IF
  100 CONTINUE

      DO 120 IVAR = IXLHS, MAXVAR
            DO 110 ITM = 1, 3
               IF (ITM .NE. ICURTM) THEN
                  IF (ISTVAR(ITM,IVAR) .GT. 0) THEN
                     IFROM = ISTVAR(ICURTM,IVAR)
                     ITO = ISTVAR(ITM,IVAR)
                     CALL CPYVAR (TYPVAR(IVAR), MAXNE,
     &                            VARVAL(1,IFROM), VARVAL(1,ITO))
                  END IF
               END IF
  110       CONTINUE
  120 CONTINUE

      RETURN
      END
