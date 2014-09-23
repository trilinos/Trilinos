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

c$$$      program test
c$$$      parameter (n=20)
c$$$      dimension a(n)
c$$$      integer indx(n)
c$$$
c$$$      integer*4 seed
c$$$      seed=0
c$$$      
c$$$      call srand(time())
c$$$      do i=1, n
c$$$        a(i) = rand(seed)
c$$$      end do
c$$$
c$$$      do i=1,n
c$$$        write (*,*) a(i)
c$$$      end do
c$$$      write (*,*) '-------'
c$$$
c$$$      call indexx(a, indx, n, .true.)
c$$$      
c$$$      do i=1,n
c$$$        write (*,*) i, indx(i), a(indx(i))
c$$$      end do
c$$$      
c$$$      end
C
C------------------------------------------------------------------------
C     SUBROUTINE INDEXN: Indexes an array ARRAY, that is
C           it outputs an array INDX such that ARRAY(INDX(J)) is in 
C           ascending order for J=1,2,...,N.  The input quantities N and 
C           ARRAY are not changed.
C
C     ARRAY (NROW, *)  -  Array to be sorted, sorted on row IROW
C     NROW             -  Row dimension of ARRAY
C     IROW             -  Row of ARRAY to be sorted
C     INDX  (modified) -  Sorted order of ARRAY
C     N                -  Number of elements in ARRAY
C     INIT             -  .FALSE. if INDX already setup
C                         .TRUE.  if INDX must be initialized
C------------------------------------------------------------------------
C
      subroutine indexn (a, nrow, irow, indx, n, init)
      
      dimension a(nrow, *)
      integer indx(0:*)
      integer n
      logical init
      
      integer start, bottom
      integer temp
 
      if (init) then
        do i=1, n
            indx(i-1) = i
          end do
      end if

      do start = (n-2)/2, 0, -1
        call siftdownn(a, nrow, irow, indx, start, n)
      end do
 
      do bottom = n-1, 1, -1
        temp = indx(0)
        indx(0) = indx(bottom)
        indx(bottom) = temp
        call siftdownn(a, nrow, irow, indx, 0, bottom)
      end do
      end
 
      subroutine siftdownn(a, nrow, irow, indx, start, bottom)
 
      real a(nrow,*)
      integer indx(0:*)
      
      integer start, bottom
      integer child, root
      integer temp
 
      root = start
      do while(root*2 + 1 .lt. bottom)
        child = root * 2 + 1
 
        if ((child + 1 .lt. bottom) .and.
     *    (a(irow, indx(child)) .lt. a(irow, indx(child+1)))) then
          child = child + 1
        end if
 
        if (a(irow, indx(root)) .lt. a(irow, indx(child))) then
          temp = indx(child)
          indx(child) = indx(root)
          indx(root) = temp
          root = child
        else
          return
        end if  
      end do    
      return
      end


