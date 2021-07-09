C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION CKCROS (IMID, IH, NLNKF, LINKF1, XN, YN, ZN)
C=======================================================================

C   --*** CKCROS *** (MESH) Hide partial lines by face
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --CKCROS deletes partial lines that have a visible node along the edge
C   --of a quarilateral and the hidden node is within the area defined
C   --by the vectors from the corner to the two adjacent nodes.  Such
C   --a face totally hides the partial line.
C   --
C   --Parameters:
C   --   IMID     - IN - the node in the face that is closest to the
C   --                   visible node
C   --   IH       - IN - the hidden nodes
C   --   NLNKF    - IN - the number of nodes per face
C   --   LINKF1   - IN - the connectivity for the face containing the edge
C   --   XN,YN,ZN - IN - the nodal coordinates

      include 'debug.blk'
C     common /debugc/ cdebug
C     common /debugn/ idebug
C     character*8 cdebug

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)

      CKCROS = .FALSE.

C   --Find left and right node for face

      DO 100 ILINK = 1, 4
         IF (LINKF1(ILINK) .EQ. IMID) IXMID = ILINK
  100 CONTINUE
      IF (IXMID .LT. 4) THEN
         ILFT = LINKF1(IXMID+1)
      ELSE
         ILFT = LINKF1(1)
      END IF
      IF (IXMID .GT. 1) THEN
         IRGT = LINKF1(IXMID-1)
      ELSE
         IRGT = LINKF1(4)
      END IF
      IF ((ILFT .EQ. IH) .OR. (IRGT .EQ. IH)) GOTO 110
      XMID = XN(IMID)
      YMID = YN(IMID)
      XLFT = XN(ILFT)
      YLFT = YN(ILFT)
      XRGT = XN(IRGT)
      YRGT = YN(IRGT)
      IF ((XLFT-XMID)*YRGT + (XRGT-XLFT)*YMID
     &   + (XMID-XRGT)*YLFT .LT. 0.0) THEN
         I = IRGT
         IRGT = ILFT
         ILFT = I
         XLFT = XN(ILFT)
         YLFT = YN(ILFT)
         XRGT = XN(IRGT)
         YRGT = YN(IRGT)
      END IF

C   --Check 2 triangles formed by 2 adjacent nodes of the face
C   --and the node have a positive area

      X0 = XN(IH)
      Y0 = YN(IH)
      ALFT = (XLFT-XMID)*Y0 + (X0-XLFT)*YMID + (XMID-X0)*YLFT
      IF (ALFT .LT. 0.0) GOTO 110
      ARGT = (XMID-XRGT)*Y0 + (X0-XMID)*YRGT + (XRGT-X0)*YMID
      IF (ARGT .LT. 0.0) GOTO 110

C   --Check if the node's Z is behind the face's Z at that point

      ZMIN = MIN
     &   (ZN(LINKF1(1)), ZN(LINKF1(2)), ZN(LINKF1(3)), ZN(LINKF1(4)))

      IF (ZMIN .LE. ZN(IH)) THEN

C      --Calculate normal at visible node with diagonal node and
C      --left or right node

         IF (IXMID .GT. 2) THEN
            IDIA = LINKF1(IXMID-2)
         ELSE
            IDIA = LINKF1(IXMID+2)
         END IF

         ZMID = ZN(IMID)

         AMIN = MIN (ALFT, ARGT)

         IF (AMIN .EQ. ALFT) THEN
            ISIDA = ILFT
            ISIDB = IDIA
         ELSE IF (AMIN .EQ. ARGT) THEN
            ISIDA = IDIA
            ISIDB = IRGT
         END IF
         AX = XN(ISIDA) - XMID
         AY = YN(ISIDA) - YMID
         AZ = ZN(ISIDA) - ZMID
         BX = XN(ISIDB) - XMID
         BY = YN(ISIDB) - YMID
         BZ = ZN(ISIDB) - ZMID
         XNORM = 0.5 * (AY*BZ - BY*AZ)
         YNORM = 0.5 * (AZ*BX - BZ*AX)
         ZNORM = 0.5 * (AX*BY - BX*AY)

         IF ((XNORM*(X0-XMID) + YNORM*(Y0-YMID)
     &      + ZNORM*(ZN(IH)-ZMID)) .GE. 0.0) GOTO 110
      END IF

C   --Mark totally hidden lines

      CKCROS = .TRUE.

  110 CONTINUE
      RETURN
      END
