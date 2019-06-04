C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Log: paint.f,v $
C Revision 1.3  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:06:40  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/08/04  15:58:41  gdsjaar
c Performed some optimization to speed up blot - Inlined the locrea function
c
c Revision 1.2  1990/12/14  08:54:31  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE PAINT (VARNP, LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, XF, YF, ZF, ISVOK, FMIN, FMAX, *)
C=======================================================================

C   --*** PAINT *** (DETOUR) Paint contours
C   --   Modified by John H. Glick - 10/26/88
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --PAINT paints contour sections in a color sequence.  The contour
C   --algorithm assumes that the elements do not contain internal nodes.
C   --The element interpolation field is approximated by logically drawing
C   --lines from each node to the element center and, thusly, dividing the
C   --element into triangles.  Contour sections are then drawn by connecting
C   --the intersection points of the sub-element edges and the contour
C   --plane.
C   --
C   --The element block status indicates which element blocks are active
C   --for contours.  No contours are drawn in inactive elements.
C   --
C   --Parameters:
C   --   VARNP - IN - the contour function values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NXFAC - IN - the number of ordered faces
C   --   IXFAC - IN - the indices of the ordered faces
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   ISVOK - IN - ISVOK(i) is true iff the contour variable is defined
C   --      for element block i (always true if nodal variable)
C   --   FMIN, FMAX - IN - the minimum and maximum function value
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses NCNTR, NOCMIN, NOCMAX of /CNTR/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM
      COMMON /CNTR/   CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC,
     &   CINTV(256), NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX
      LOGICAL CINTOK, LINCON, NOCMIN, NOCMAX

      REAL VARNP(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XF(*), YF(*), ZF(*)
      LOGICAL ISVOK(NELBLK)

      LOGICAL GRABRT

C   --Fill the CINTV array with the contour intervals

      IF (.NOT. CINTOK) THEN
         DO 100 NC = 1, NCNTR+1
            CINTV(NC) = CNTRI (NC)
  100    CONTINUE
      END IF

C   --Set the contour minimum and maximum (in case NOCMIN, etc.)
      CNTRMN = FMIN - ABS (DELC)
      CNTRMX = FMAX + ABS (DELC)

      DO 120 IX = 1, NXFAC
         IFAC = IXFAC(IX)
         IELB = 0
         IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
         IF (NLNKF(IELB) .LE. 2) GOTO 120

         IF ((.NOT. IS3DIM)  .AND.  (NLNKF(IELB) .EQ. 9)) THEN
            NNPF = 8
         ELSE
            NNPF = NLNKF(IELB)
         ENDIF

         IF (ISVOK(IELB)) THEN

C         --Compute the minimum and maximum values for the face
            CALL CFMAX (VARNP, NLNKF(IELB), LINKF(IXL), FEMIN, FEMAX)

C         --Find limiting contour intervals
              difx = abs(cintv(1) - femax)
              difn = abs(cintv(1) - femin)
              imaxc = 1
              iminc = 1
              do 10 i=2, ncntr+1
                difxx = abs(cintv(i) - femax)
                if (difx .gt. difxx) then
                  difx = difxx
                  imaxc = i
                end if
                difnn = abs(cintv(i) - femin)
                if (difn .gt. difnn) then
                  difn = difnn
                  iminc = i
                end if
 10           continue
            IF (DELC .GE. 0.0) THEN
c$$$               MINC = LOCREA (FEMIN, NCNTR+1, CINTV)
c$$$               MAXC = LOCREA (FEMAX, NCNTR+1, CINTV)
              minc = iminc
              maxc = imaxc
              IF (FEMIN .LT. CINTV(MINC)) MINC = MINC - 1
              IF (FEMAX .LT. CINTV(MAXC)) MAXC = MAXC - 1
            ELSE
c$$$               MINC = LOCREA (FEMAX, NCNTR+1, CINTV)
c$$$               MAXC = LOCREA (FEMIN, NCNTR+1, CINTV)
              minc = imaxc
              maxc = iminc
              IF (FEMAX .GE. CINTV(MINC)) MINC = MINC - 1
              IF (FEMIN .GE. CINTV(MAXC)) MAXC = MAXC - 1
            END IF

            IF (NOCMIN .AND. (MINC .LT. 1)) THEN
               MINC = 1
               IF (MINC .GT. MAXC) MAXC = MINC
            END IF
            IF (NOCMAX .AND. (MAXC .GT. NCNTR)) THEN
               MAXC = NCNTR
               IF (MINC .GT. MAXC) MINC = MAXC
            END IF
            IF ((MAXC .LT. 1) .OR. (MINC .GT. NCNTR)) THEN
               MINC = -1
               MAXC = -1
            END IF

         ELSE
C         --Not selected, paint face black
            MINC = -1
            MAXC = -1
         END IF

C      --Skip this contour if the values are outside the contour range

         IF (MINC .EQ. MAXC) THEN

C         --Face is entirely inside one contour area

            CALL GRCOLR (MINC)
            CALL SOLIDF (NNPF, LINKF(IXL), XN, YN, ZN)

         ELSE

C         --Face has several contour areas

            IF ((MINC .LT. 1) .OR. (MAXC .GT. NCNTR)) THEN

C            --Some part of face is outside contour limits, paint black

               CALL GRCOLR (-1)
               MINC = MAX (MINC, 1)
               MAXC = MIN (MAXC, NCNTR)
            END IF

            DO 110 NC = MINC, MAXC
               CALL GRCOLR (NC)
               IF (DELC .GE. 0.0) THEN
                  CNTR0 = CNTRI (NC)
                  IF (NOCMIN .AND. (NC .EQ. 1)) CNTR0 = CNTRMN
                  CNTR1 = CNTRI (NC+1)
                  IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR1 = CNTRMX
               ELSE
                  CNTR1 = CNTRI (NC)
                  IF (NOCMIN .AND. (NC .EQ. 1)) CNTR1 = CNTRMX
                  CNTR0 = CNTRI (NC+1)
                  IF (NOCMAX .AND. (NC .EQ. NCNTR)) CNTR0 = CNTRMN
               END IF

               IF ((FEMIN .GE. CNTR1) .OR. (FEMAX .LT. CNTR0))
     &            print *, femin, cntr1, femax, cntr0

               IF (GRABRT ()) RETURN 1
               CALL PAINTF (CNTR0, CNTR1,
     &            VARNP, NLNKF(IELB), LINKF(IXL),
     &            XN, YN, ZN, XF(IFAC), YF(IFAC), ZF(IFAC))
  110       CONTINUE
         END IF
  120 CONTINUE

      CALL PLTFLU

      RETURN
      END
