C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
C=======================================================================
      SUBROUTINE INVCON (ia, NELBLK, IDELB, NUMELB, NUMLNK, LINK, NUMNP,
     *  ICON, ICSCR, NODIND, lisnp, nout, mapn, mape, domapn, domape,
     $     DOBLK, DOELE, ebtype)
C=======================================================================

C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements in all blocks
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   LINK - IN - the connectivity array

      include 'params.blk'
      include 'dbase.blk'

      INTEGER IA(*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      INTEGER ICON(NELBLK,*)
      INTEGER ICSCR(2,NELBLK)
      INTEGER NODIND(2,NUMNP)
      INTEGER LISNP(0:*)
      INTEGER MAPN(*)
      INTEGER MAPE(*)
      LOGICAL DOMAPN
      LOGICAL DOMAPE
      LOGICAL FIRST, DOBLK, DOELE
      CHARACTER*(MXSTLN) EBTYPE(*)
      DATA FIRST /.TRUE./
      SAVE index, ndcon
      
      if (first) then
        first = .false.
        call iniint (numnp*nelblk, 0, icon)
        call iniint (2*numnp, 0, nodind)

        ISLNK = 1
        DO IELB = 1, NELBLK
          if (ebtype(ielb) .eq. 'nsided' .or.
     *      ebtype(ielb) .eq. 'NSIDED') THEN
           CALL invcn1 (IELB, 1, NUMLNK(IELB), 
     &          LINK(ISLNK), icon, nelblk)
           ISLNK = ISLNK + NUMLNK(IELB) 
          ELSE
           CALL invcn1 (IELB, NUMELB(IELB), NUMLNK(IELB), 
     &          LINK(ISLNK), icon, nelblk)
           ISLNK = ISLNK + NUMLNK(IELB) * NUMELB(IELB)
         END IF
        end do

C     ... Build up a node index map so we can store the inverse node
C     -element connectivity in a linear array.
        index = 1
        do i = 1, numnp
           isum = 0
           do j=1, nelblk
              isum = isum + icon(j,i)
           end do
           nodind(1,i) = index
           nodind(2,i) = 0
           index = index + isum
        end do
        call mdrsrv('NODCON', ndcon, index)
        ISLNK = 1
        IELST = 0
        DO IELB = 1, NELBLK
          if (ebtype(ielb) .eq. 'nsided' .or.
     *      ebtype(ielb) .eq. 'NSIDED') THEN
            call mdrsrv('NNPE', nnpe, numelb(ielb))
            CALL EXGECPP(NDB, EXEBLK, IDELB(ielb), ia(nnpe), IERR)
            call invcn2n(IELST, numelb(ielb), ia(nnpe), 
     &        LINK(ISLNK), nodind, numnp, ia(ndcon))
            call mddel('NNPE')
            ISLNK = ISLNK + NUMLNK(IELB)
            IELST = IELST + NUMELB(IELB)
         ELSE
           call invcn2(IELST, NUMELB(IELB), NUMLNK(IELB), 
     &          LINK(ISLNK), nodind, numnp, ia(ndcon))
           ISLNK = ISLNK + NUMLNK(IELB) * NUMELB(IELB)
           IELST = IELST + NUMELB(IELB)
         END IF
        end do
        if (domape) then
           do i = 0, index-1
              ia(ndcon+i) = mape(ia(ndcon+i))
           end do
        end if
      end if
      
C=======================================================================
C ... OUTPUT
C=======================================================================
      
      if (nout .gt. 0) then
        if (domapn) write (nout, 10005) 'Node'
      else
        if (domapn) write (*, 10005) 'Node'
      end if

      if (nout .gt. 0) then
        if (domape) write (nout, 10005) 'Element'
      else
        if (domape) write (*, 10005) 'Element'
      end if

C ... Node/Block inverse connectivity      
      if (doblk) then
         if (nout .gt. 0) then
            WRITE (nout, 10000)
         else
            WRITE (*, 10000)
         end if            
         do i = 1, lisnp(0)
            INP = LISNP(I)
            
C     ... Get nonzero blocks for this node.
            icnt = 0
            do j=1, nelblk
               if (icon(j,inp) .gt. 0) then
                  icnt = icnt + 1
                  icscr(1,icnt) = idelb(j)
                  icscr(2,icnt) = icon(j,inp)
               end if
            end do
            
            if (domapn) then
               id = mapn(inp)
            else
               id = inp
            end if
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10020) ID, (icscr(1,j), icscr(2,j),j=1,icnt)
            ELSE
               WRITE (*, 10020)    ID, (icscr(1,j), icscr(2,j),j=1,icnt)
            END IF
         end do
      end if

C ... Node/Element inverse connectivity      
      if (doele) then
         if (nout .gt. 0) then
            WRITE (nout, 10010)
         else
            WRITE (*, 10010)
         end if            
         do i = 1, lisnp(0)
            INP = LISNP(I)
            
            ibeg = ndcon + nodind(1,inp) - 2
            icnt = nodind(2,inp)
            if (domapn) then
               id = mapn(inp)
            else
               id = inp
            end if
            if (domape) then
              IF (NOUT .GT. 0) THEN
                WRITE (NOUT, 10030) ID, (mape(ia(ibeg+k)), k=1,icnt)
              ELSE
                WRITE (*, 10030)    ID, (mape(ia(ibeg+k)), k=1,icnt)
              END IF
            else
              IF (NOUT .GT. 0) THEN
                WRITE (NOUT, 10030) ID, (ia(ibeg+k), k=1,icnt)
              ELSE
                WRITE (*, 10030)    ID, (ia(ibeg+k), k=1,icnt)
              END IF
           end if
         end do
      end if

10000 FORMAT (/,1x, 'Inverse Node-Block Connectivity:',/,
     *  ' Output is: Block ID:# times in this block')
10010 FORMAT (/,1x, 'Inverse Node-Element Connectivity:',/,
     *  ' Output is: id of elements connected to node')
10005 FORMAT (1X, A,' ids are Global')
10020 FORMAT (1X, 'Node', I10,2X, 5 (1X, i10,':',i2.2), :, /,
     &  (17X, 5(1X, i10,':',i2.2)))
10030 FORMAT (1X, 'Node', I10,':',2X, 8 (1X, i10), :, /,
     &  (17X, 8(1X, i10)))
      RETURN
      END
      
      subroutine invcn1(ielb, numelb, numlnk, link, icon, nelblk)

      integer link(numlnk,*)
      integer icon(nelblk,*)
      
      do i=1, numelb
         do j = 1, numlnk
            node = link(j,i)
            icon(ielb,node) = icon(ielb,node) + 1
         end do
      end do
      return
      end
      
      subroutine invcn2(ielst, numelb, numlnk, link, 
     $     nodind, numnp, nodcon)
      integer link(numlnk,*)
      integer nodind(2,numnp)
      integer nodcon(*)
      
      do i=1, numelb
        do j = 1, numlnk
          node = link(j,i)
          index = nodind(1,node)
          icnt  = nodind(2,node)
          nodcon(index+icnt) = ielst+i
          nodind(2,node) = nodind(2,node) + 1
       end do
      end do
      return
      end

      subroutine invcn2n(ielst, numelb, nnpe, link, 
     $     nodind, numnp, nodcon)
      integer nnpe(*), link(*)
      integer nodind(2,numnp)
      integer nodcon(*)
      
      ind = 0
      do i=1, numelb
        numlnk = nnpe(i)
        do j = 1, numlnk
          node = link(ind + j)
          index = nodind(1,node)
          icnt  = nodind(2,node)
          nodcon(index+icnt) = ielst+i
          nodind(2,node) = nodind(2,node) + 1
        end do
        ind = ind + numlnk
      end do
      return
      end



