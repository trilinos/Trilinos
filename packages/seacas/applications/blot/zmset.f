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

C $Log: zmset.f,v $
C Revision 1.3  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1994/07/21 15:28:22  gdsjaar
C Moved more commons into includes.
C
c Revision 1.1  1994/04/07  20:18:28  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.2  1990/12/14  08:59:46  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ZMSET (XZMMIN, XZMMAX, YZMMIN, YZMMAX, XN, YN,
     &   LINSET, IPSET, NPART)
C=======================================================================

C   --*** ZMSET *** (MESH) Delete partial lines outside zoom window
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --ZMSET deletes partial lines with both nodes outside the zoom window,
C   --and the line not crossing the zoom window.
C   --
C   --Parameters:
C   --   XZMMIN, XZMMAX, YZMMIN, YZMMAX - IN - the box enclosing the
C   --      zoom window
C   --   XN, YN - IN - the nodal coordinates
C   --   LINSET - IN/OUT - the sorted line set
C   --   IPSET - IN/OUT - the indices of the partial line set
C   --   NPART - IN/OUT - the number of lines in the partial line set
C   --
C   --Common Variables:
C   --   Uses LLNSET of /D3NUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'd3nums.blk'

      REAL XN(*), YN(*)
      INTEGER LINSET(LLNSET,*)
      INTEGER IPSET(*)

      YLN = YZMMAX - YZMMIN

      nhid = 0
      IP = 1
  100 CONTINUE
      IF (IP .LE. NPART) THEN

C      --Delete partial line if both nodes are outside zoom window

         N1 = LINSET(1,IPSET(IP))
         X1 = XN(N1)
         Y1 = YN(N1)
         IF ((X1 .LE. XZMMIN) .OR. (X1 .GE. XZMMAX)
     &      .OR. (Y1 .LE. YZMMIN) .OR. (Y1 .GE. YZMMAX)) THEN
            N2 = LINSET(2,IPSET(IP))
            X2 = XN(N2)
            Y2 = YN(N2)
            IF ((X2 .LE. XZMMIN) .OR. (X2 .GE. XZMMAX)
     &         .OR. (Y2 .LE. YZMMIN) .OR. (Y2 .GE. YZMMAX)) THEN
               IF (((X1 .LE. XZMMIN) .AND. (X2 .LE. XZMMIN)) .OR.
     &            ((X1 .GE. XZMMAX) .AND. (X2 .GE. XZMMAX)) .OR.
     &            ((Y1 .LE. YZMMIN) .AND. (Y2 .LE. YZMMIN)) .OR.
     &            ((Y1 .GE. YZMMAX) .AND. (Y2 .GE. YZMMAX))) THEN
                  IPSET(IP) = IPSET(NPART)
                  NPART = NPART - 1
                  IP = IP - 1
                  nhid = nhid + 1

               ELSE
C               --Calculate the intersection of the window and the partial line
C               --Solve the simultaneous equations:
C               --   X = X2 + (X1 - X2) * TVH = XLIM
C               --   Y = Y2 + (Y1 - Y2) * TVH = YZMMIN + (YZMMAX - YZMMIN) * TLN

                  IF (X1 .LE. XZMMIN) THEN
                     XLIM = XZMMIN
                  ELSE
                     XLIM = XZMMAX
                  END IF

                  XVH = X1 - X2
                  YVH = Y1 - Y2
                  XLH = XLIM - X2
                  YLH = YZMMIN - Y2
                  IF ((XVH .NE. 0.0) .AND. (YLN .NE. 0.0)) THEN
                     TLN = - (-YVH * XLH + XVH * YLH) / (XVH * YLN)
                     IF ((TLN .LE. 0) .OR. (TLN .GE. 1)) THEN
                        IPSET(IP) = IPSET(NPART)
                        NPART = NPART - 1
                        IP = IP - 1
                        nhid = nhid + 1
                     END IF
                  END IF
               END IF
            END IF
         END IF

         IP = IP + 1
         GOTO 100
      END IF
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'partial lines outside zoom =', nhid

      RETURN
      END
