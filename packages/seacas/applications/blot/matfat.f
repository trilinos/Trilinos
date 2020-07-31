C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION MATFAT (LINKF1, MAXNPF, NPFS, iel, numnp, IERR)
C=======================================================================

C   --*** MATFAT *** (MESH) Match face with existing faces
C   --   Written by Amy Gilkey - revised 02/19/88
C   --   Revised by John Glick - 10/20/88
C   --              Sam Key, 06/01/85
C   --

C .... MODIFIED FOR TRIANGULAR FACES OF TET ELEMENTS....

C   --MATFAT searches for the given face in the list of previously
C   --identified faces.  If a match is found, a pointer to the
C   --matching face is returned.
C   --
C   --The search is done by looking at the list of the faces which
C   --contain each node (the NPFS array).  The lists of two diagonally
C   --opposite nodes in the given face are searched for a reference
C   --to a common face, which is a "match".  The references to the
C   --common face are deleted for all nodes in the face.
C   --
C   --Parameters:
C   --   LINKF1 - IN - the nodes of the face
C   --   MAXNPF - IN - the maximum length of the NPFS entry
C   --   NPFS - IN/OUT - the list of unmatched faces containing a node;
C   --      (i,0) = the length of the list
C   --   IEL - IN - the element containing the face (for error message)
C   --   IERR - OUT - = 0 if routine executes with no errors
C   --                = 1 if an error was detected.

      include 'minmax.blk'
      INTEGER LINKF1(4)
      INTEGER NPFS(NUMNP, 0:MAXNPF)

      MATFAT = 0
      IERR = 0

      IF ((NPFS(LINKF1(1),0) .GT. 0)
     &   .AND. (NPFS(LINKF1(2),0) .GT. 0)
     &   .AND. (NPFS(LINKF1(3),0) .GT. 0)) THEN

C      --Check all prior faces using this node, by looking for the
C      --diagonally opposite node of this face

         INF1 = LINKF1(1)
         INF3 = LINKF1(3)
         DO 150 I1 = 1, NPFS(INF1,0)
            IOLDF = NPFS(INF1,I1)
            DO 140 I3 = 1, NPFS(INF3,0)
               IF (IOLDF .EQ. NPFS(INF3,I3)) THEN

C               --First two nodes matche so check other node

                  INF2 = LINKF1(2)
                  DO 100 I2 = 1, NPFS(INF2,0)
                     IF (IOLDF .EQ. NPFS(INF2,I2)) GOTO 110
  100             CONTINUE
                  GOTO 140

  110             CONTINUE

C               --Matching faces so point to face and zero out
C               --nodal references to face

                  L = NPFS(INF1,0)
                  IF (L .GT. I1) NPFS(INF1,I1) = NPFS(INF1,L)
                  NPFS(INF1,0) = L - 1

                  L = NPFS(INF3,0)
                  IF (L .GT. I3) NPFS(INF3,I3) = NPFS(INF3,L)
                  NPFS(INF3,0) = L - 1

                  L = NPFS(INF2,0)
                  IF (L .GT. I2) NPFS(INF2,I2) = NPFS(INF2,L)
                  NPFS(INF2,0) = L - 1

                  MATFAT = IOLDF
                  maxnod = max(maxnod, inf1, inf2, inf3)
                  minnod = min(minnod, inf1, inf2, inf3)

                  GO TO 160
               END IF
  140       CONTINUE
  150    CONTINUE
      END IF

  160 CONTINUE
      RETURN
      END
