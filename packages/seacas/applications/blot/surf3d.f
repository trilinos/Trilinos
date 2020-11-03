C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SURF3D (LENSC, NLNKSC, LINKSC, IF2ESC,
     &   LENF, NLNKF, LINKF, IF2EL, IF2EL2, NFACES, LENLNK)
C=======================================================================

C   --*** SURF3D *** (MESH) Group 3D faces
C   --   Written by Amy Gilkey - revised 03/29/88
C   --
C   --SURF3D groups all the faces in the 3D mesh into one of the
C   --following parts:
C   --   n) Surface faces of element block n
C   --   NELBLK+1) Interior faces
C   --
C   --Parameters:
C   --   LENSC - IN - the cumulative face counts by element block for LINKSC
C   --   NLNKSC - IN - the number of nodes per face for LINKSC
C   --   LINKSC - IN - the unsorted connectivity for all faces
C   --   IF2ESC - IN - the element number of each face in LINKSC;
C   --      IF2ESC(2,x) = 0 iff surface face
C   --      IF2ESC(2,x) > 0 iff interior face
C   --      IF2ESC(2,x) < 0 iff duplicate interior face
C   --   LENF - OUT - the cumulative face counts by element block
C   --   NLNKF - OUT - the number of nodes per face
C   --   LINKF - OUT - the connectivity for all faces
C   --   IF2EL - OUT - the element number of each face
C   --   IF2EL2 - OUT - the secondary element number of each face
C   --   NFACES - OUT - the number of returned faces
C   --   LENLNK - OUT - the length of the returned LINKF array
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENSC(0:NELBLK)
      INTEGER NLNKSC(NELBLK)
      INTEGER LINKSC(*)
      INTEGER IF2ESC(NDIM-1,*)
      INTEGER LENF(0:NELBLK+4)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)

C   --Move surface faces into appropriate element block

      LENF(0) = 0
      NFACES = 0
      IXL0 = 0
      DO 110 IELB = 1, NELBLK
         NLNKF(IELB) = NLNKSC(IELB)
         DO 100 IFAC = LENSC(IELB-1)+1, LENSC(IELB)
            IF (IF2ESC(2,IFAC) .EQ. 0) THEN
               IXS = IDBLNK (IELB, IFAC, LENSC, NLNKSC)
               NFACES = NFACES + 1
               CALL CPYINT (NLNKF(IELB), LINKSC(IXS), LINKF(IXL0+1))
               IF2EL(NFACES) = IF2ESC(1,IFAC)
               IF2EL2(NFACES) = 0
               IXL0 = IXL0 + NLNKF(IELB)
            END IF
  100    CONTINUE
         LENF(IELB) = NFACES
  110 CONTINUE

C   --Move interior faces

      DO 130 IELB = 1, NELBLK
         DO 120 IFAC = LENSC(IELB-1)+1, LENSC(IELB)
            IF (IF2ESC(2,IFAC) .GT. 0) THEN
               NFACES = NFACES + 1
               IXS = IDBLNK (IELB, IFAC, LENSC, NLNKSC)
               CALL CPYINT (NLNKSC(IELB), LINKSC(IXS), LINKF(IXL0+1))
               IF2EL(NFACES) = IF2ESC(1,IFAC)
               IF2EL2(NFACES) = IF2ESC(2,IFAC)
               IXL0 = IXL0 + NLNKF(IELB)
            END IF
  120    CONTINUE
  130 CONTINUE
      LENF(NELBLK+1) = NFACES

      LENF(NELBLK+2) = NFACES
      LENF(NELBLK+3) = NFACES
      LENF(NELBLK+4) = NFACES

      LENLNK = IXL0
      RETURN
      END
