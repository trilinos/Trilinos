*-----------------------------------------------------------------------
*\Documentation
*
*\Name: TSPLOOPS
*
*\Description:
*
*  Driver routine which reads in a set of Harwell-Boeing format
*  matrices and test a set of sparse loops
*
*\Usage:
*     tsploops < in
*
*\Input files:
*    Unit 5   Reads the following information from standard input:
*
*             1) NSETS, the number of data sets to be processed.
*                These must be stored in Harwell-Boeing format.
*
*             2) The names the NSETS data files to be
*                read in by SCSCRD.  Each name should be on a separate
*                line.
*
*\Remarks:
*     None.
*
*\Examples:
*     None.
*
*\Enddoc
*
*---------------------------------------------------------------------
*
*\Lib
*
*\Local variables:
*     xxxxxx  integer
*
*\Author:
*     Mike Heroux
*     Cray Research, Inc.
*
*\References:
*     None.
*
*\Keywords:
*     sparse BLAS; CSC, CSR
*
*\Routines called:
*
*\Revision history:
*     08/09/92: Original Implementation (Mike Heroux)
*
*\Implementation Details:
*     None.
*
*\Performance Data:
*     None.
*
*\Endlib
*
*---------------------------------------------------------------------
      program tsploops
*
*     -----------------------------
*     Specifications for parameters
*     -----------------------------
      implicit none
      integer nnzmax, mmax, nrhsmx
      parameter ( 
     &    nnzmax   = 3 700 000,   ! Max number of entries in matrix
     &    mmax     =    47 000,   ! Max problem dimension
     &    nrhsmx   =      5)   ! Max problem dimension
      logical debug
      parameter(debug=.false.)
c
      integer indx (nnzmax),  pntr  (mmax+1)
      integer indxt(nnzmax),  pntrt (mmax+1)
      integer rhsptr(mmax+1), rhsind(mmax)
      integer ntrials, nsets, iset, mtxunt, job, m, n, nnz, nrhs
      integer ierr, i, j, ibgn, iend, itrials
      real*8  fdiv, t1, tmmc1, tmmc2, tmmr1, tmmr2
      real*8  fmflops_mmc1, fmflops_mmc2, fmflops_mmr1, fmflops_mmr2
      real*8  err_mm1, err_mm2, errall, fnops
      real*8  val   (nnzmax), valt  (nnzmax)
      real*8  x(mmax*nrhsmx), rhs (mmax*nrhsmx)
      real*8 soln(mmax*nrhsmx)
      real*8 rhsc1(mmax*nrhsmx), rhsc2(mmax*nrhsmx)
      real*8 rhsr1(mmax*nrhsmx), rhsr2(mmax*nrhsmx)
c
      character*80 title
      character*80 filname
      character*8  key
      character*3  type, rhstyp
      real*8 norm, second
      external norm, second
*
*     --------------------------
*     First executable statement
*     --------------------------
c
c     ------------
c     Print header
c     ------------
      write(*,*)
      write(*,100)
     &'Description',      'Sparse Column',         'Sparse Row'
      write(*,110)
     &'Data Set ','   Dim   ',' # entries  ',
     &'V1 MFLP','V2 MFLP',
     &'V1 MFLP','V2 MFLP','Max Err ','NRHS '
      write(*,110)
     &'---------','---------','------------',
     &'-------','-------',
     &'-------','-------','--------',' ----'

 100  format(' ',t11,a11,t35,a15,t52,a15)
 110  format(' ',a8,t11,a9 ,t21,a12,
     &t36,a7,t44,a7,
     &t54,a7,t62,a7,t72,a8,t80,a5)
c
c.....Read the following from Fortran unit 5:
c
c     1) NSETS, the number of data sets to be processed.
c        These must be stored in Harwell-Boeing format.
c
c     2) The names the NSETS data files to be
c        read in by SCSCRD.  Each name should be on a separate
c        line.
c
c.....get NSETS
      read (5,*)nsets
c
c.....Now loop through the data sets
      do 9001 iset = 1, nsets
c
c.....get filename
      read (5,'(a)') filname
c
c.....open file
      mtxunt = 1
      open (unit=mtxunt, file=filname, status='old', err=9000)
      if (debug) write(*,*)'INFO:  Opened file: ',filname
c
c     ----------------------------
c     Read Harwell-Boeing data set
c     ----------------------------
c.....read matrix structure only
      job = 1
      call scscrd (mmax, nnzmax, job, mtxunt, val, indx, pntr,
     &             rhs, rhsptr, rhsind, nrhs, soln, x, 
     &             rhstyp, m, n, nnz, title, 
     &             key, type, ierr)
      if (ierr .ne. 0) then
         print *, ' Error from SCSCRD: ', ierr, ' for ', key
         go to 9000
      end if
c
c.....close file
      close(mtxunt)
      if (debug) write(*,*)'INFO:  Closed file: ',filname
        
c
c.....Make up fake values (diag = 1).
      if (job.lt.2) then
         if (debug) write(*,*)'INFO:  Generating values and solution'
         do 5 j = 1, n
            ibgn = pntr(j)
            iend = pntr(j+1) - 1
            val(ibgn) = 1.0D0
            ibgn = ibgn + 1
            if (ibgn.le.iend) then
               fdiv = - 0.5D0/float(iend-ibgn+1)/float(2*n)
               do 10 i = ibgn, iend
                  val(i) = fdiv * (float(j) + float(indx(i)))
 10            continue
            endif
 5       continue
      endif
c
c.....initialize timer
      t1 = second()

c
c     -------------------
c     Form transpose of L
c     -------------------
      call csrcsc(n, 1, 1, val, indx, pntr, valt, indxt, pntrt )
      if (debug) write(*,*)'INFO:  Converted to sparse row format'

      do 101 nrhs = 1, nrhsmx
c
c     -----------------------
c     define solution
c     -----------------------
        do 15 i=1,n*nrhs
          soln(i) = float(i)
C           soln(i) = 1.0
 15     continue

      ntrials = max(30,min(100000000/(n*nrhs),100))
      fnops = 2.0 * float(nrhs)  * float(nnz) * float(ntrials)
c
c     *********************
c     Sparse Column Version 1 
c     *********************
      t1 = second()
      do itrials = 1,ntrials
      call scscmm1(n, n, val, indx, pntr, soln, rhsc1, nrhs)
      end do
      tmmc1 = second() - t1
      fmflops_mmc1 = fnops/tmmc1*1.0D-6
      if (debug) write(*,*)'INFO: Sparse column MM interleaved'
c
c     *********************
c     Sparse Column Version 2 
c     *********************
      t1 = second()
      do itrials = 1,ntrials
      call scscmm2(n, n, val, indx, pntr, soln, rhsc2, nrhs)
      end do
      tmmc2 = second() - t1
      fmflops_mmc2 = fnops/tmmc2*1.0D-6
      if (debug) write(*,*)'INFO: Sparse column MM non-interleaved'
c
c     *********************
c     Sparse Row Version 1 
c     *********************
      t1 = second()
      do itrials = 1,ntrials
      call scsrmm1(n, n, valt, indxt, pntrt, soln, rhsr1, nrhs)
      end do
      tmmr1 = second() - t1
      fmflops_mmr1 = fnops/tmmr1*1.0D-6
      if (debug) write(*,*)'INFO:  Sparse Row MM interleaved'
c
c     *********************
c     Sparse Row Version 2 
c     *********************
      t1 = second()
      do itrials = 1,ntrials
      call scsrmm2(n, n, valt, indxt, pntrt, soln, rhsr2, nrhs)
      end do
      tmmr2 = second() - t1
      fmflops_mmr2 = fnops/tmmr2*1.0D-6
      if (debug) write(*,*)'INFO:  Sparse Row MM non-interleaved'
c
c.....Verify solution version 1
      call swaxpy(n*nrhs, -1.0D0, rhsc1, rhsr1, rhs)
      err_mm1 = norm( n*nrhs, rhs )
      if (debug) write(*,*)'INFO: Error Col/Row V1', err_mm1
c
c.....Verify solution version 2
      call swaxpy(n*nrhs, -1.0D0, rhsc2, rhsr2, rhs)
      err_mm2 = norm( n*nrhs, rhs )
      if (debug) write(*,*)'INFO: Error Col/Row V2', err_mm2
      errall = max(err_mm1, err_mm2)
c
c     -------------------------------
c     Print results for this data set
c     -------------------------------
      write(*,120)key,n, nnz,
     &fmflops_mmc1,fmflops_mmc2,
     &fmflops_mmr1,fmflops_mmr2,errall,nrhs

 120  format(' ',a8,t11,i9 ,t21,i12,
     &t36,f7.2,t44,f7.2,
     &t54,f7.2,t62,f7.2,t72,E8.2,t80,i5)
 101  continue
      goto 9001
c
 9000 continue
      write(*,*)' Error reading file'
 9001 continue
      end
