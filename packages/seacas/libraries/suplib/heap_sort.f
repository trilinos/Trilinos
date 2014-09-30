C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
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

      subroutine heapsort(a,n)
      
      dimension a(0:*)
      integer n
      
      integer start, bottom
      real temp
 
      do start = (n-2)/2, 0, -1
        call hs_siftdown(a, start, n)
      end do
 
      do bottom = n-1, 1, -1
        temp = a(0)
        a(0) = a(bottom)
        a(bottom) = temp
        call hs_siftdown(a, 0, bottom)
      end do
      end
 
      subroutine hs_siftdown(a, start, bottom)
 
      real a(0:*)
      integer start, bottom
      integer child, root
      real temp
 
      root = start
      do while(root*2 + 1 .lt. bottom)
        child = root * 2 + 1
 
        if ((child + 1 .lt. bottom) .and.
     *    (a(child) .lt. a(child+1))) then
          child = child + 1
        end if
 
        if (a(root) .lt. a(child)) then
          temp = a(child)
          a(child) = a (root)
          a(root) = temp
          root = child
        else
          return
        end if  
      end do    
      return
      end


