C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION NZAWAY (NLNKF, LINKF1, XN, YN, ZN, HIDENP)
C=======================================================================

C   --*** NZAWAY *** (MESH) Determine if warped face is hidden
C   --   Written by Amy Gilkey - revised 05/02/88
C   --
C   --NZAWAY determines how many nodes of a possible warped face
C   --on a 3D surface are hidden.  A node is hidden if and only if the
C   --outward normal formed by the node and the two adjacent nodes points
C   --"into" the plotting surface.  The nodes are marked as visible if
C   --they point away from the plotting surface.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates; the Z-coordinate is
C   --      pointing towards the viewer (out of the plotting plane)
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      INTEGER HIDENP(*)

      NZAWAY = 0
      IMID = LINKF1(3)
      IRGT = LINKF1(4)
      AX = XN(IRGT) - XN(IMID)
      AY = YN(IRGT) - YN(IMID)
      IF ((AX .NE. 0.0) .OR. (AY .NE. 0.0)) THEN
         ISAME = 0
      ELSE
         IMID = LINKF1(2)
         IRGT = LINKF1(3)
         BX = XN(IRGT) - XN(IMID)
         BY = YN(IRGT) - YN(IMID)
         IMID = LINKF1(3)
         IRGT = LINKF1(4)
         ISAME = IMID
      END IF

      DO 100 ILINK = 1, 4
         IMID = IRGT
         IRGT = LINKF1(ILINK)

C      --Form X and Y components of corner of face

         IF (ISAME .EQ. 0) THEN
            BX = -AX
            BY = -AY
         END IF
         AX = XN(IRGT) - XN(IMID)
         AY = YN(IRGT) - YN(IMID)

C      --Form Z component of normal vector to corner, and make node
C      --visible if normal points forward

         IF ((AX .NE. 0.0) .OR. (AY .NE. 0.0)) THEN
            IF (AX*BY .LE. BX*AY) THEN
               NZAWAY = NZAWAY + 1
               IF (ISAME .NE. 0) NZAWAY = NZAWAY + 1
            ELSE
               IXRGT = ILINK
               HIDENP(IMID) = KNVIS
               IF (ISAME .NE. 0) HIDENP(ISAME) = KNVIS
            END IF
            ISAME = 0
         ELSE
            ISAME = IMID
         END IF
  100 CONTINUE

      IF (NZAWAY .EQ. 1) THEN

C      --If single node visible, make two adjoining nodes visible

         IF (IXRGT .GT. 2) THEN
            IXLFT = IXRGT - 2
         ELSE
            IXLFT = IXRGT + 2
         END IF
         HIDENP(LINKF1(IXLFT)) = KNVIS
         HIDENP(LINKF1(IXRGT)) = KNVIS
      END IF

      if ((nzaway .ne. 4) .and. (nzaway .ne. 0)) then
         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 2))
     &      write (*, '(1x,a,i5,2x,4i5)') 'warped', nzaway, linkf1
      end if

      RETURN
      END
