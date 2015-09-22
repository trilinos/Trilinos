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

C $Log: zmedge.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:18:24  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:59:43  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ZMEDGE (XZMMIN, XZMMAX, YZMMIN, YZMMAX, XN, YN,
     &   IEDSET, NEDGES)
C=======================================================================

C   --*** ZMEDGE *** (MESH) Delete edges outside zoom window
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --ZMEDGE deletes edges with both nodes outside the zoom window,
C   --and the line not crossing the zoom window.
C   --
C   --Parameters:
C   --   XZMMIN, XZMMAX, YZMMIN, YZMMAX - IN - the box enclosing the
C   --      zoom window
C   --   XN, YN - IN - the nodal coordinates
C   --   IEDSET - IN/OUT - the edge line set;
C   --      (0) = face defining edge; 0 to delete edge
C   --   NEDGES - IN - the number of lines in the edge set

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      REAL XN(*), YN(*)
      INTEGER IEDSET(0:2,*)

      nhid = 0
      DO 100 IEDG = 1, NEDGES
         IF (IEDSET(0,IEDG) .EQ. 0) GOTO 100

C      --Delete edge if both nodes are outside zoom window

         N1 = IEDSET(1,IEDG)
         X1 = XN(N1)
         Y1 = YN(N1)
         IF ((X1 .LE. XZMMIN) .OR. (X1 .GE. XZMMAX)
     &      .OR. (Y1 .LE. YZMMIN) .OR. (Y1 .GE. YZMMAX)) THEN
            N2 = IEDSET(2,IEDG)
            X2 = XN(N2)
            Y2 = YN(N2)
            IF (((X1 .LE. XZMMIN) .AND. (X2 .LE. XZMMIN)) .OR.
     &         ((X1 .GE. XZMMAX) .AND. (X2 .GE. XZMMAX)) .OR.
     &         ((Y1 .LE. YZMMIN) .AND. (Y2 .LE. YZMMIN)) .OR.
     &         ((Y1 .GE. YZMMAX) .AND. (Y2 .GE. YZMMAX))) THEN
               IEDSET(0,IEDG) = 0
               nhid = nhid + 1
            END IF
         END IF

  100 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'edges outside zoom window =', nhid

      RETURN
      END
