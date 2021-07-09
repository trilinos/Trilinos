C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TRNXYZ (XN, YN, XN3, YN3, ZN3, IXNP, NRNP, ZCORD)
C=======================================================================

C   --*** NEWXYZ *** (GEN3D) Calculate 3D coordinates
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --
C   --NEWXYZ calculates the coordinate array for the 3D database.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --   NPCEN - IN - the node numbers of the center nodes by column and row
C   --   ZCORD - SCRATCH - size = NNREPL, holds z coordinate for transformations
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      INTEGER IXNP(*), NRNP(*)
      REAL ZCORD(NNREPL)

C   --Copy Y coordinate from original

         DO 520 INP = 1, NUMNP
            JNP0 = IXNP(INP) - 1
            DO 510 NR = 1, NRNP(INP)
               YN3(JNP0+NR) = YN(INP)
  510       CONTINUE
  520    CONTINUE

C      --For translations, repeat the X coordinate and add a Z coordinate

      IBLK = 0
      NXTNR = 1

      ZEND = 0.0
    1    CONTINUE
         IBLK = IBLK + 1
         IF (NRTRAN(IBLK) .GT. 0) THEN
            ZBEG = ZEND
            ZEND = ZBEG + D3TRAN(IBLK)
            CALL INIGRD (-ZBEG, -ZEND, ZGRAD(IBLK),
     *         NRTRAN(IBLK), NRTRAN(IBLK)+1, ZCORD(NXTNR) )
            NXTNR = NXTNR + NRTRAN(IBLK)
            IF (IBLK .LT. MAXINT) GO TO 1
         END IF

C      --Repeat the X coordinate and add the calculated Z coordinate

      DO 30 INP = 1, NUMNP
         JNP = IXNP(INP) - 1
         DO 20 NR = 1, NNREPL
            XN3(JNP+NR) = XN(INP)
            ZN3(JNP+NR) = ZCORD(NR)
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
