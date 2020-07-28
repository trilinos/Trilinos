C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TRNXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
C=======================================================================

C   --*** TRNXYZ *** (GENSHELL) Calculate 3D coordinates for translation
C   --
C   --TRNXYZ calculates the coordinate array for the 3D database.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
C   --   Uses XOFFS, YOFFS, ZOFFS of /XYZOFF/
C   --   Uses ROT3D, ROTMAT of /XYZROT/

      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(NUMEL)

C   --Copy X and Y coordinates from original, set Z equal to 0.0

      DO 10 INP = 1, NUMNP
         XN3(INP) = XN(INP)
         YN3(INP) = YN(INP)
         ZN3(INP) = 0.0
 10   CONTINUE

C   --Set the element attributes equal to the thickness

      DO 20 IEL = 1, NUMEL
         ATRIB(IEL) = DIM3
 20   CONTINUE

      RETURN
      END
