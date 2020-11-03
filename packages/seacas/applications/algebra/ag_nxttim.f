C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NXTTIM ()
C=======================================================================

C   --*** NXTTIM *** (ALGEBRA) Move current time variables to last time
C   --   Written by Amy Gilkey - revised 12/10/87
C   --
C   --NXTTIM moves the values of the variables for the current time to the
C   --location for the last time.  The move is implemented by swapping
C   --the location pointers for the last and current time.
C   --
C   --Parameters:
C   --
C   --Common Variables:
C   --   Sets ISTVAR of /VAR../
C   --   Uses NUMINP, IXLHS of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      include 'ag_namlen.blk'
      include 'ag_var.blk'

c      LOGICAL WHOTIM

c      LOGICAL ISHIST

C   --If previous time requested, swap pointers for current and previous time

      DO 100 IVAR = 1, NUMINP
c         ISHIST = ((TYPVAR(IVAR) .EQ. 'H') .OR. (TYPVAR(IVAR) .EQ. 'T'))
c         IF (ISHIST .OR. WHOTIM) THEN
         IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
            I = ISTVAR(ICURTM,IVAR)
            ISTVAR(ICURTM,IVAR) = ISTVAR(ILSTTM,IVAR)
            ISTVAR(ILSTTM,IVAR) = I
         END IF
c         END IF
  100 CONTINUE

      DO 110 IVAR = IXLHS, MAXVAR
c         ISHIST = ((TYPVAR(IVAR) .EQ. 'H') .OR. (TYPVAR(IVAR) .EQ. 'T'))
c         IF (ISHIST .OR. WHOTIM) THEN
         IF (ISTVAR(ILSTTM,IVAR) .GT. 0) THEN
            I = ISTVAR(ICURTM,IVAR)
            ISTVAR(ICURTM,IVAR) = ISTVAR(ILSTTM,IVAR)
            ISTVAR(ILSTTM,IVAR) = I
         END IF
c         END IF
  110 CONTINUE

      RETURN
      END
