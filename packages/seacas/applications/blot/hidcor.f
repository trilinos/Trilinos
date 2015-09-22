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

C=======================================================================
      SUBROUTINE HIDCOR (LINSET, IPSET, NPART, IEDSET, NEDGES,
     &                   LENF, NLNKF, LINKF, XN, YN, ZN,
     *                   MREF, LREF)
C=======================================================================

C   --*** HIDCOR *** (MESH) Hide partial lines by face
C   --   Written by Amy Gilkey - revised 10/22/87
C   --
C   --HIDCOR deletes partial lines that have a visible node as a corner
C   --of a quarilateral and the hidden node is within the area defined
C   --by the vectors from the corner to the two adjacent nodes.  Such
C   --a face totally hides the partial line.
C   --
C   --Parameters:
C   --   LINSET   - IN  - the sorted line set
C   --   IPSET    - I/O - the indices of the partial line set
C   --   NPART    - I/O - the number of lines in the partial line set
C   --   IEDSET   - I/O - the edge line set;
C   --                    (0) = face defining edge; 0 to delete edge
C   --   NEDGES   - I/O - the number of lines in the edge set
C   --   LENF     - IN  - the cumulative face counts by element block
C   --   NLNKF    - IN  - the number of nodes per face
C   --   LINKF    - IN  - the connectivity for all faces
C   --   XN,YN,ZN - IN  - the nodal coordinates
C   --
C   --Common Variables:
C   --   Uses LLNSET of /D3NUMS/

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LINSET(LLNSET,*)
      INTEGER IPSET(*)
      INTEGER IEDSET(0:2,*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER MREF(*), LREF(*)
      
      LOGICAL CKCROS

C   --Check each edge line against each partial line for overlap

      NOLDPT = NPART

      call iniint(numnpf, 0, mref)
      call iniint(numnpf, 0, lref)

      do ip=1,npart
        lref(LINSET(1,IPSET(IP))) = ip
        lref(LINSET(2,IPSET(IP))) = ip
      end do

      do ip=npart,1,-1
        mref(LINSET(1,IPSET(IP))) = ip
        mref(LINSET(2,IPSET(IP))) = ip
      end do

      nhid = 0
      DO 120 IEDG = 1, NEDGES
         IFAC = IEDSET(0,IEDG)
         N1 = IEDSET(1,IEDG)
         N2 = IEDSET(2,IEDG)

         imin = max(1, mref(n1), mref(n2))
         imax = min(npart, lref(n1), lref(n2))
         do ip = imin, imax
           if (ipset(ip) .gt. 0) then
C         --Process only if visible node is corner node of face
             IV = LINSET(1,IPSET(IP))
             IF ((N1 .EQ. IV) .OR. (N2 .EQ. IV)) THEN
               IH = LINSET(2,IPSET(IP))
               IF (N2 .EQ. IH) GOTO 110
               
C            --Check if face hides the visible part of the line
               IELB = 0
               IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
               IF (CKCROS (IV, IH, NLNKF(IELB), LINKF(IXL),
     &           XN, YN, ZN)) THEN
                 
C               --Move totally hidden lines from partial set to totally
C               --hidden lines
                 ipset(ip) = -ipset(ip)
                 nhid = nhid + 1
                 
               END IF
             END IF
           end if
 110       CONTINUE
         end do
 120  CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &  write (*, '(1x,a,i5)') 'partial lines hidden by face =', nhid

      nhid = 0
      ip = 1
 100  continue
      if (ip .le. npart) then
        if (ipset(ip) .lt. 0) then
          I = -IPSET(IP)
          IPSET(IP) = IPSET(NPART)
          IPSET(NPART) = I
          NPART = NPART - 1
          IP = IP - 1
          nhid = nhid + 1
        end if
        ip = ip + 1
        go to 100
      end if
C   --Delete the edges which are totally hidden lines

      nhid = 0
      DO 140 IP = NPART+1, NOLDPT
         IV = LINSET(1,IPSET(IP))
         IH = LINSET(2,IPSET(IP))
         DO 130 IEDG = 1, NEDGES
            IF (IEDSET(1,IEDG) .EQ. IV) THEN
               IF (IEDSET(2,IEDG) .EQ. IH) THEN
                  IEDSET(0,IEDG) = 0
                  nhid = nhid + 1
                  GOTO 140
               END IF
            END IF
  130    CONTINUE
  140 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'edges hidden by face =', nhid

      RETURN
      END
