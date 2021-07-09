C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNPLOT (NPTIMS, NPTIMW, XLN, YLN, ZLN, IXNODE, *)
C=======================================================================

C   --*** LNPLOT *** (PATHLN) Plot the pathlines
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNPLOT plots the pathlines.  For 3D meshes, the pathlines are rotated.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of points on a history pathline
C   --   NPTIMW - IN - the number of points on a non-history pathline
C   --   XLN, YLN, ZLN - IN - the pathline data
C   --   IXNODE - SCRATCH - size = NPTIMS (3D only)
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NLNCRV, ILVID of /LNVARS/
C   --   Uses ROTMAT, ROTCEN of /ROTOPT/

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'lnvars.blk'
      include 'd3nums.blk'
      include 'rotopt.blk'

      REAL XLN(NPTIMS,NLNCRV), YLN(NPTIMS,NLNCRV), ZLN(NPTIMS,NLNCRV)
      INTEGER IXNODE(NPTIMS)

      LOGICAL GRABRT
      LOGICAL NUMCRV
      CHARACTER TYP
      CHARACTER*8 LABSID

      lintyp = 1
      isytyp = 0
      numcrv = .false.
      labsid = 'NONE'

      IF (IS3DIM) THEN
         DO 100 I = 1, NPTIMS
            IXNODE(I) = I
  100    CONTINUE
      END IF

      DO 110 NP = 1, NLNCRV

         IF (GRABRT()) RETURN 1
         CALL GRCOLR (NP)
         CALL GRSYMB (LINTYP, ISYTYP, NP)

         CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
         IF (TYP .EQ. 'H') THEN
            NPTS = NPTIMS
         ELSE
            NPTS = NPTIMW
         END IF

C      --Bl_Rotate the 3D pathline

         IF (IS3DIM) THEN
            CALL BL_ROTATE (NPTS, IXNODE, ROTMAT, ROTCEN,
     &         XLN, YLN, ZLN, XLN, YLN, ZLN)
         END IF

C      --Plot pathline

         IF (GRABRT()) RETURN 1
         CALL PLTCUR (XLN(1,NP), YLN(1,NP), NPTS)

         IF (NUMCRV) THEN
            IF (GRABRT()) RETURN 1
            CALL GRNCRV (LABSID, NP, NPTS,
     &         XLN(1,NP), YLN(1,NP), (LINTYP .EQ. 0))
         END IF

C      --Finish plot

         CALL PLTFLU
  110 CONTINUE

      RETURN
      END
