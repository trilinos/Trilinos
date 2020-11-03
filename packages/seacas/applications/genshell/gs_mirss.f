C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MIRSS (IDFRO, IDBCK, NLINK,
     &   NSSUR, NSSFRO, NSSBCK, LTSES3)
C=======================================================================

C   --*** MIRSS *** (GEN3D) Modifies sideset node order to account
C                           for mirroring about axes
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from WRESS written by Amy Gilkey
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface side sets; (0) = length
C   --   IDBCK - IN - ids for back surface side sets; (0) = length
C   --   NSSUR - IN - the number of nodes in the surface side set
C   --   NSSFRO - IN - the nodes in the front surface side set
C   --   NSSBCK - IN - the nodes in the back surface side set
C   --   LTSES3 - IN - the sides for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NUMESS of /DBNUMS/
C   --   Uses LESSNO of /DBNUM3/
C   --
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER NSSFRO(NLINK,*), NSSBCK(NLINK,*)
      INTEGER LTSES3(*)

      LOGICAL ANYESS

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYESS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMESS .GT. 0)

C   --Write 3D

      IF (ANYESS) THEN

c$$$        IF (NLINK .EQ. 4) THEN
c$$$          INCR = 2
c$$$        ELSE IF (NLINK .EQ. 8 .OR. NLINK .EQ. 9) THEN
c$$$          INCR = 3
c$$$        END IF
c$$$C   --Sidesets on 'side' of shell mesh (original 2d sidesets) - lines 2 nodes
c$$$        DO 10 NL = 1, LESSNO, INCR
c$$$          LNTMP = LTNES3(NL)
c$$$          LTNES3(NL) = LTNES3(NL+1)
c$$$          LTNES3(NL+1) = LNTMP
c$$$ 10     CONTINUE

C   --Sidesets on front surface of shell mesh - faces 4 nodes
        IF (MOD(NSSUR, NLINK) .NE. 0) THEN
          CALL PRTERR('FATAL',
     *      'Incorrect front surface node list in mirss')
          STOP 'MIRSS'
        END IF
        NFACE = NSSUR / NLINK
        DO 20 NL = 1, NFACE
          LNTMP = NSSFRO(2, NFACE)
          NSSFRO(2, NFACE) = NSSFRO(4, NFACE)
          NSSFRO(4, NFACE) = LNTMP
          LNTMP = NSSBCK(2, NFACE)
          NSSBCK(2, NFACE) = NSSBCK(4, NFACE)
          NSSBCK(4, NFACE) = LNTMP
 20     CONTINUE

        IF (NLINK .EQ. 8 .OR. NLINK .EQ. 9) THEN
          DO 30 NL = 1, NFACE
            LNTMP = NSSFRO(5, NFACE)
            NSSFRO(5, NFACE) = NSSFRO(8, NFACE)
            NSSFRO(8, NFACE) = LNTMP
            LNTMP = NSSFRO(6, NFACE)
            NSSFRO(6, NFACE) = NSSFRO(7, NFACE)
            NSSFRO(7, NFACE) = LNTMP
            LNTMP = NSSBCK(5, NFACE)
            NSSBCK(5, NFACE) = NSSBCK(8, NFACE)
            NSSBCK(8, NFACE) = LNTMP
            LNTMP = NSSBCK(6, NFACE)
            NSSBCK(6, NFACE) = NSSBCK(7, NFACE)
            NSSBCK(7, NFACE) = LNTMP
 30       CONTINUE
        END IF
      END IF

      RETURN
      END

