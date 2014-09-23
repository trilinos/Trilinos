C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

      subroutine getssn(ia, ierr)

      include 'dbase.blk'

      integer ia(*)
      
      CALL MDFIND ('LTSESS', KLTSSS, LESSEL)
      CALL MDFIND ('IXNESS', KIXNSS, NUMESS)
      CALL MDFIND ('IDESS',  KIDSS,  NUMESS)
      CALL MDFIND ('LTNNSS', KLTNNN, LESSEL)
      CALL MDFIND ('LTNESS', KLTNSS, LESSNL)
      CALL MDFIND ('NNESS',  KNNSS,  NUMESS)

c ...Convert sides to nodes.... a(kltsss), 
C offset into element list for current side set
      isoff = 0
C     node count for current side set
      nodcnt = 0
      do i=0,numess-1
C     update index array            
         ia(kixnss+i)=nodcnt+1
C     get num of sides & df            
         call exgsp(ndb,ia(kidss+i),nsess,ndess,ierr)
      if (ierr .gt. 0) goto 170 
         
C     get side set nodes
         call exgssn(ndb,ia(kidss+i),ia(kltnnn+isoff),
     &        ia(kltnss+nodcnt),ierr) 
      if (ierr .gt. 0) goto 170
      nness = 0
C     sum node counts to calculate next index
         do ii=0,nsess-1 
            nness=nness+ia(kltnnn+isoff+ii)
         end do
         ia(knnss+i)=nness
         nodcnt=nodcnt+nness
         isoff=isoff+nsess
      end do
      
 170  continue
      return
      end
      
