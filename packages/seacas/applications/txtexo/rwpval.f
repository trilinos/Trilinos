C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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




