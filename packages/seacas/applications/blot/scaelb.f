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

C $Log: scaelb.f,v $
C Revision 1.3  2009/03/25 12:36:47  gdsjaar
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
C Revision 1.1  1994/04/07 20:10:33  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:48  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCAELB (A, USESEL, IELBST,
     &   VALMN, NUMMN, XYZMN, ISTMN, VALMX, NUMMX, XYZMX, ISTMX,
     &   VALMIN, NUMMIN, XYZMIN, ISTMIN, VALMAX, NUMMAX, XYZMAX, ISTMAX)
C=======================================================================

C   --*** SCAELB *** (BLOT) Scale element variable for selected blocks
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAELB returns the minimum and maximum values for the selected
C   --element blocks.  The minimum and maximums for each block have already
C   --been calculated and stored.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   USESEL - IN - use the element blocks selected array iff true,
C   --      else all selected
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   VALMN, VALMX - IN - the minimum and maximum value for each
C   --      element block
C   --   NUMMN, NUMMX - IN - the element number of the minimum and
C   --      maximum value for each element block
C   --   XYZMN, XYZMX - IN - the coordinates of NUMMN, NUMMX
C   --   ISTMN, ISTMX - IN - the step number of the minimum and maximum
C   --      value for each element block
C   --   VALMIN, VALMAX - OUT - the minimum and maximum value
C   --      (with selected element block and birth/death)
C   --   NUMMIN, NUMMAX - OUT - the element number of the minimum and
C   --      maximum value
C   --   XYZMIN, XYZMAX - OUT - the coordinates of NUMMIN, NUMMAX
C   --   ISTMIN, ISTMAX - OUT - the step number of the minimum and maximum value
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK, NVAREL, NSTEPS of /DBNUMS/

      include 'dbnums.blk'

      DIMENSION A(*)
      LOGICAL USESEL
      INTEGER IELBST(*)
      REAL VALMN(0:NELBLK), VALMX(0:NELBLK)
      INTEGER NUMMN(0:NELBLK), NUMMX(0:NELBLK)
      REAL XYZMN(3,0:NELBLK), XYZMX(3,0:NELBLK)
      INTEGER ISTMN(0:NELBLK), ISTMX(0:NELBLK)
      REAL XYZMIN(3), XYZMAX(3)

      LOGICAL INIT

C   --Calculate min/max for element by element block

      IF (USESEL) THEN
         INIT = .TRUE.
         DO 100 IELB = 1, NELBLK
            IF ((.NOT. USESEL) .OR. (IELBST(IELB) .GT. 0)) THEN
               IF (ISTMN(IELB) .GT. 0) THEN
                  IF (INIT) THEN
                     IMIN = IELB
                     IMAX = IELB
                     INIT = .FALSE.
                  ELSE
                     IF (VALMN(IMIN) .GT. VALMN(IELB)) THEN
                        IMIN = IELB
                     END IF
                     IF (VALMX(IMAX) .LT. VALMX(IELB)) THEN
                        IMAX = IELB
                     END IF
                  END IF
               END IF
            END IF
  100    CONTINUE
      ELSE
         INIT = .FALSE.
         IMIN = 0
         IMAX = 0
      END IF

      IF (INIT) THEN
         VALMIN = 0.0
         NUMMIN = 0
         DO 110 I = 1, 3
            XYZMIN(I) = 0.0
  110    CONTINUE
         ISTMIN = 0
         VALMAX = 0.0
         NUMMAX = 0
         DO 120 I = 1, 3
            XYZMAX(I) = 0.0
  120    CONTINUE
         ISTMAX = 0
      ELSE
         VALMIN = VALMN(IMIN)
         NUMMIN = NUMMN(IMIN)
         DO 130 I = 1, 3
            XYZMIN(I) = XYZMN(I,IMIN)
  130    CONTINUE
         ISTMIN = ISTMN(IMIN)
         VALMAX = VALMX(IMAX)
         NUMMAX = NUMMX(IMAX)
         DO 140 I = 1, 3
            XYZMAX(I) = XYZMX(I,IMAX)
  140    CONTINUE
         ISTMAX = ISTMX(IMAX)
      END IF

      RETURN
      END
