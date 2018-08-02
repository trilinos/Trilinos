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

C $Log: conlab.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:57:10  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:58  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CONLAB (ICNTR, CNTR, NHIT, LABINC, LINSET,
     &   VARNP, XN, YN, ZN, IN2ELB, *)
C=======================================================================

C   --*** CONLAB *** (DETOUR) Label contour line
C   --   Written by Amy Gilkey - revised 06/01/87
C   --   D. P. Flanagan, 03/30/83
C   --
C   --CONLAB labels a contour line that intersects the given line with the
C   --contour letter.  If LABINC is greater than one, only every LABINCth
C   --crossing is labeled.
C   --
C   --Parameters:
C   --   ICNTR - IN - the contour label number
C   --   CNTR - IN - the contour value
C   --   NHIT - IN/OUT - the number of "hits" so far, initialize to zero
C   --   LABINC - IN - the number of "hits" to skip before labelling
C   --   LINSET - IN - the line set
C   --   VARNP - IN - the contour function values
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (DXO = .0015, DYO = .0015)

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LINSET(3)
      REAL VARNP(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IN2ELB(*)

      LOGICAL INTERP_BL
      LOGICAL GRABRT
      LOGICAL EXISTS

      EXISTS (M) = (MOD(M,2) .NE. 0)

      IF (IS3DIM) THEN
C      --Skip invisible line
         IF (LINSET(3) .EQ. 0) GOTO 100
      END IF

      N1 = LINSET(1)
      N2 = LINSET(2)
      IF ((IN2ELB(N1) .GE. 0) .AND. (IN2ELB(N2) .GE. 0)) THEN

         IF (INTERP_BL (CNTR, VARNP(N1), VARNP(N2), PSI)) THEN
            NHIT = NHIT + 1
            IF (NHIT .GE. LABINC) THEN
               IF (GRABRT ()) RETURN 1
               IF (IS3DIM) THEN
C               --Replace invisible node with partial line node
                  IF (ABS (LINSET(3)) .GT. 1)
     &               N2 = ABS (LINSET(3))
               END IF
               X0 = XN(N1) * (1.-PSI) + XN(N2) * PSI
               Y0 = YN(N1) * (1.-PSI) + YN(N2) * PSI
               CALL MP2PT (1, X0, Y0, DX0, DY0, MASK)
               IF (EXISTS (MASK)) CALL PLTXTS (DX0+DXO, DY0+DYO,
     &            CHAR (ICHAR('A')-1+ICNTR))
               NHIT = 0
            END IF
         END IF
      END IF

  100 CONTINUE
      RETURN
      END
