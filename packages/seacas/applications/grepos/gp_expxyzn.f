C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EXPXYZN (XN, YN, ZN, XEXPL, YEXPL, ZEXPL,
     &   NUMNPS, IDNPS, NNNPS, IXNNPS, LTNNPS, NUMNP, NDIM,
     $     MODE)
C=======================================================================

C   --*** EXPXYZ *** (GREPOS) Modify coordinates for each block separately
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --EXPXYZ modifies the coordinate array for the database.
C   --       each block is treated separately if not connected
C   --
C   --Parameters:
C   --   XN, YN, ZN - OUT - the coordinates
C   --   MODE  - 1 = Explode
C   --           2 = Scale
C   --           3 = Randomize
C   --           4 = Node Randomize

      REAL XN(*), YN(*), ZN(*)
      REAL XEXPL(*), YEXPL(*), ZEXPL(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      INTEGER MODE

      if (MODE .NE. 4) RETURN

C ... Randomize Each Nodeset node
      IDUM = 1
      DO 150 I = 1, NUMNPS
        if (XEXPL(I) .ne. 0.0 .or. YEXPL(I) .ne. 0.0 .or.
     *    ZEXPL(I) .ne. 0.0) THEN
          DO 140 INOD = IXNNPS(I), IXNNPS(I)+nnnps(i)-1
            NODE = LTNNPS(INOD)
            XN(NODE) = (2.0*RAN1(IDUM)-1.0) * XEXPL(I) + XN(NODE)
            YN(NODE) = (2.0*RAN1(IDUM)-1.0) * YEXPL(I) + YN(NODE)
            IF (NDIM .EQ. 3) THEN
               ZN(NODE) = (2.0*RAN1(IDUM)-1.0) * ZEXPL(I) + ZN(NODE)
             END IF
 140       CONTINUE
         END IF
 150   CONTINUE
      RETURN
      END
