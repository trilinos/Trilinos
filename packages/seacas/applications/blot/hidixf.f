C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDIXF (HIDEF, HIDENP, LENF, NLNKE, NLNKF, LINKF,
     &   ZF, NXFAC, IXFAC, IIXFAC, NAMELB)
C=======================================================================

C   --*** HIDIXF *** (MESH) Order visible faces by Z coordinate
C   --   Written by Amy Gilkey - revised 02/29/88
C   --
C   --HIDIXF orders the visible faces by the Z coordinate of the face center.
C   --The ordering is expressed by indexing the faces from minimum Z to
C   --maximum Z.
C   --
C   --Parameters:
C   --   HIDEF - IN - face status (as in HIDDEN)
C   --   HIDENP - IN - node status (as in HIDDEN)
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   ZF - IN - the face Z coordinates
C   --   NXFAC - OUT - the number of ordered faces
C   --   IXFAC - OUT - the indices of the ordered faces
C   --   IIXFAC - SCRATCH - size = number of faces
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      PARAMETER (KFVIS=0, KFNODH=10, KFPOUT=20, KFOUT=90, KFAWAY=100)
      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER HIDEF(*)
      INTEGER HIDENP(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKE(NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      REAL ZF(*)
      INTEGER IXFAC(*)
      INTEGER IIXFAC(*)
      CHARACTER*(*) NAMELB(*)

      NXFAC = 0
      NNXFAC = 0
      DO 120 IELB = 1, NELBLK
         IF (NLNKE(IELB) .GT. 4 .OR. NAMELB(IELB)(:3) .EQ. 'TET') THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 110 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LT. KFAWAY) THEN
                  NHID = 0
                  DO 100 ILINK = 1, NLNKF(IELB)
                     IF (HIDENP(LINKF(IXL0+ILINK)) .GT. KNVIS)
     &                  NHID = NHID + 1
  100             CONTINUE
                  IF (NHID .LE. 0) THEN
                     NNXFAC = NNXFAC + 1
                     IIXFAC(NNXFAC) = IFAC
                  ELSE
                     NXFAC = NXFAC + 1
                     IXFAC(NXFAC) = IFAC
                  END IF
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
  110       CONTINUE
         END IF
         IF (NLNKE(IELB) .EQ. 4 .AND. NAMELB(IELB)(:3) .NE. 'TET') THEN
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 1110 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LT. KFAWAY) THEN
                  NHID = 0
                  DO 1100 ILINK = 1, NLNKF(IELB)
                     IF (HIDENP(LINKF(IXL0+ILINK)) .GT. KNVIS)
     &                  NHID = NHID + 1
 1100             CONTINUE
                  IF (NHID .LE. 0) THEN
                     NNXFAC = NNXFAC + 1
                     IIXFAC(NNXFAC) = IFAC
                  ELSE
                     NXFAC = NXFAC + 1
                     IXFAC(NXFAC) = IFAC
                  END IF
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
 1110       CONTINUE
         END IF
  120 CONTINUE

      call indexx(zf, ixfac, nxfac, .FALSE.)
c$$$      DO 140 IX = 1, NXFAC
c$$$         IMIN = IX
c$$$         ZMIN = ZF(IXFAC(IMIN))
c$$$         DO 130 IX1 = IX+1, NXFAC
c$$$            IF (ZF(IXFAC(IX1)) .LT. ZMIN) THEN
c$$$               IMIN = IX1
c$$$               ZMIN = ZF(IXFAC(IMIN))
c$$$            END IF
c$$$  130    CONTINUE
c$$$         ISAV = IXFAC(IMIN)
c$$$         IXFAC(IMIN) = IXFAC(IX)
c$$$         IXFAC(IX) = ISAV
c$$$  140 CONTINUE

      DO 150 IX = 1, NNXFAC
         NXFAC = NXFAC + 1
         IXFAC(NXFAC) = IIXFAC(IX)
  150 CONTINUE

      RETURN
      END
