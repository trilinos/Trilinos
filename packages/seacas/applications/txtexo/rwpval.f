C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 

      subroutine rwpval(ntxt, ndb, a, ia, c, nelblk, numnps, numess,
     &  idelb, idnps, idess, *)
      include 'exodusII.inc'
      real a(*)
      integer ia(*)
      character*1 c(*)
      integer idelb(*), idnps(*), idess(*)
      integer cerr
      
C ... Skip comment record
      read (ntxt, *, end=100, err=100)
      
C ... Element block properties
      read (ntxt, *, end=110, err=110) numebp

      call mcrsrv ('EBPNAM', iebpn, numebp * mxstln)
      call mdrsrv ('EBPVAL', iebpv, numebp * nelblk)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         return 1
       end if
       
       call rwpval1(ntxt, ndb, EXEBLK, numebp, nelblk, idelb,
     &   c(iebpn), a(iebpv), *200)
       
       call mcdel ('EBPNAM')
       call mddel ('EBPVAL')
      call mdstat (nerr, mem)
      if (nerr .ne. 0) then
         call memerr()
         return 1
      end if

C ... Node set properties
      read (ntxt, *, end=120, err=120) numnsp
      call mcrsrv ('NSPNAM', inspn, numnsp * mxstln)
      call mdrsrv ('NSPVAL', inspv, numnsp * numnps)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         return 1
      end if

      call rwpval1 (ntxt, ndb, EXNSET, numnsp, numnps, idnps, 
     &             c(inspn), a(inspv), *200)

      call mcdel ('NSPNAM')
      call mddel ('NSPVAL')
      call mdstat (nerr, mem)
      if (nerr .ne. 0) then
         call memerr()
         return 1
      end if
      
C ... Side set properties      
      read (ntxt, *, end=120, err=120) numssp
      call mcrsrv ('SSPNAM', isspn, numssp * mxstln)
      call mdrsrv ('SSPVAL', isspv, numssp * numess)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         return 1
      end if

      call rwpval1 (ntxt, ndb, EXSSET, numssp, numess, idess,
     &             c(isspn), a(isspv), *200)

      call mcdel ('SSPNAM')
      call mddel ('SSPVAL')
      call mdstat (nerr, mem)
      if (nerr .ne. 0) then
         call memerr()
         return 1
      end if
      return
 100  continue
 110  continue
 120  continue
 200  continue
      return 1
      end



      
