C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: fixcut.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:00:59  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:40  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FIXCUT (CUTPT, CUTNRM, X, Y, Z,
     &   LENF, NLNKF, LINKF, IF2EL, IF2EL2, IE2ELB,
     &   IFACUT, IELCUT, CLASS, NEWELB)
C=======================================================================

C   --*** FIXCUT *** (MESH) Cut 3D mesh
C   --   Written by Amy Gilkey - revised 03/04/88
C   --   Revised by Ray J. Meyers, 29 May, 1990
C              modified input from three points plus a logical indicating
C              whether to reverse the implied plane normal (12 values),
C              to input of point on plane and normal of plane (6 values)
C   --
C   --FIXCUT cuts the 3D mesh along a given plane.  All faces cut by the
C   --plane or out of the plane on a cut element become surface faces.
C   --Faces that are out of the cut are moved to LENF(NELBLK+3) set.
C   --
C   --Parameters:
C   --   CUTPT - IN - a point on the cutting plane
C   --   CUTNRM - IN - the normal of the cutting plane
C   --   X, Y, Z - IN - the original nodal coordinates
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   IFACUT - SCRATCH - the face status for all faces
C   --   IELCUT - SCRATCH - the element status
C   --   CLASS - SCRATCH - integer array: for each point, the point is
C   --                     classified in, out, or on with respect to the
C                          cutting plane
C   --   NEWELB - OUT - size = LENF(NELBLK+3)
C   --
C   --Common Variables:
C   --   Uses NUMNP, NELBLK of /DBNUMS/

      PARAMETER (ISIN = +1, ISOUT = -1, ISCUT = 0, ISON = -2, ISSURF=2)

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      REAL CUTPT(3), CUTNRM(3)
      REAL X(NUMNP), Y(NUMNP), Z(NUMNP)
      INTEGER LENF(0:NELBLK+3)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER IFACUT(*), IELCUT(NUMEL)
      INTEGER CLASS(NUMNP)
      INTEGER NEWELB(*)
      LOGICAL ON, IN, OUT

C   --Rotate the coordinates to find the z-coordinates

C      DO 100 INP = 1, NUMNP
C         ZC(INP) = X(INP)*CUTMAT(1,3) + Y(INP)*CUTMAT(2,3)
C     &      + Z(INP)*CUTMAT(3,3)
C  100 CONTINUE

C
C CLASSIFY EACH POINT AS BEING ISIN, ISOUT, OR ISON
C
      DO 100 I = 1, NUMNP
         CALL CLASPT( X(I), Y(I), Z(I), CUTPT, CUTNRM, CLASS(I))
100   CONTINUE


C   --Initialize the element flags

      CALL INIINT (NUMEL, -999, IELCUT)

      DO 140 IELB = 1, NELBLK+3
         IF (IELB .LE. NELBLK) NL = NLNKF(IELB)
         IF (IELB .LE. NELBLK) IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1

         DO 130 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IELB .GT. NELBLK) NL = NLNKF(IE2ELB(IF2EL(IFAC)))

C         --Find out if the surface is in, out, on or cut

            IN = .FALSE.
            OUT = .FALSE.
            ON = .FALSE.
            DO 110 K = 1, NL
               ISTAT = CLASS(LINKF(IXL0+K))
               IF(ISTAT .EQ. ISIN) THEN
                    IN = .TRUE.
               ELSE IF(ISTAT .EQ. ISOUT) THEN
                    OUT = .TRUE.
               ELSE
                    ON = .TRUE.
               END IF
110         CONTINUE

            IF(IN) THEN
               IF(OUT) THEN
                  IFACUT(IFAC) = ISCUT
               ELSE
                  IFACUT(IFAC) = ISIN
               END IF
            ELSE IF(OUT) THEN
               IFACUT(IFAC) = ISOUT
            ELSE
               IFACUT(IFAC) = ISON
            END IF



C            ZMIN = ZC(LINKF(IXL0+1))
C            ZMAX = ZC(LINKF(IXL0+1))
C            DO 110 K = 2, NL
C               ZMIN = MIN (ZMIN, ZC(LINKF(IXL0+K)))
C               ZMAX = MAX (ZMAX, ZC(LINKF(IXL0+K)))
C  110       CONTINUE
C            IF ((ZMAX .EQ. ZCUT) .AND. (ZMIN .EQ. ZMAX)) THEN
C               IFACUT(IFAC) = ISON
C            ELSE IF (ZMAX .LT. ZCUT) THEN
C               IFACUT(IFAC) = ISIN
C            ELSE IF (ZMIN .GT. ZCUT) THEN
C               IFACUT(IFAC) = ISOUT
C            ELSE
C               IFACUT(IFAC) = ISCUT
C            END IF

C         --Set the element flag to in, out, or cut

            DO 120 K = 1, 2
               IF (K .EQ. 1) THEN
                  IEL = IF2EL(IFAC)
               ELSE
                  IEL = IF2EL2(IFAC)
               END IF
               IF (IEL .GT. 0) THEN
                  IF (IELCUT(IEL) .NE. IFACUT(IFAC)) THEN
                     IF ((IELCUT(IEL) .LE. -999)
     &                  .OR. (IFACUT(IFAC) .EQ. ISCUT)) THEN
                        IELCUT(IEL) = IFACUT(IFAC)
                     ELSE IF (((IELCUT(IEL) .EQ. ISIN)
     &                  .AND. (IFACUT(IFAC) .EQ. ISON))
     &                  .OR. ((IELCUT(IEL) .EQ. ISON)
     &                  .AND. (IFACUT(IFAC) .EQ. ISIN))) THEN
                        IELCUT(IEL) = ISCUT
                     END IF
                  END IF
               END IF
  120       CONTINUE
            IXL0 = IXL0 + NL
  130    CONTINUE
  140 CONTINUE

      DO 160 IELB = 1, NELBLK+3
         DO 150 IFAC = LENF(IELB-1)+1, LENF(IELB)

            IF (IFACUT(IFAC) .EQ. ISOUT) THEN

C            --Change an OUT face to a SURFACE face if it is part of a
C            --CUT element
               IF (IELCUT(IF2EL(IFAC)) .EQ. ISCUT) IFACUT(IFAC) = ISSURF
               IEL = IF2EL2(IFAC)
               IF (IEL .GT. 0) THEN
                  IF (IELCUT(IEL) .EQ. ISCUT) IFACUT(IFAC) = ISSURF
               END IF

            ELSE IF (IFACUT(IFAC) .EQ. ISON) THEN

C            --Change an ON face to a SURFACE face
               IFACUT(IFAC) = ISSURF

            ELSE IF (IFACUT(IFAC) .EQ. ISCUT) THEN

C            --Change a CUT face to an IN face
               IFACUT(IFAC) = ISIN

            END IF

C         --Determine the type of the new surface

            IF (IFACUT(IFAC) .EQ. ISIN) THEN

               IF (IF2EL2(IFAC) .LE. 0) THEN
C               --If surface IN face, change to a surface face
                  NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
               ELSE
C               --If interior IN face, change to an interior face
                  NEWELB(IFAC) = NELBLK+1
               END IF

            ELSE IF (IFACUT(IFAC) .EQ. ISOUT) THEN

C            --If OUT face, change to an OUT face
               NEWELB(IFAC) = NELBLK+3

            ELSE IF (IFACUT(IFAC) .EQ. ISSURF) THEN

C            --If SURFACE surface face, change to a surface face

               IF (IF2EL2(IFAC) .LE. 0) THEN
                  NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
               ELSE IF (IELCUT(IF2EL2(IFAC)) .EQ. ISOUT) THEN
                  NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
               ELSE
C               --Swap nodes to simulate surface being defined by facing element
                  NEWELB(IFAC) = - IE2ELB(IF2EL2(IFAC))
               END IF
            END IF

  150    CONTINUE
  160 CONTINUE

      RETURN
      END
