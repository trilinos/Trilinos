      program test1
      implicit none
      integer nmax, kmax
      parameter ( nmax = 250 000, kmax = 16)

      real*8 zero, one, onemill
      parameter ( zero = 0.0E0, one = 1.0E0, onemill = 1.0E6)
      real*8 c(nmax*kmax), a(nmax*kmax), b(kmax*kmax)
      real*8 c1(nmax*kmax), b1(kmax*kmax)
      character*1 trans, notrans

      integer i, n, m, k, n0, n1, k0, k1, itry, ntry
      real*8 alpha, beta, fnops, fnops150, t1, t2
      real*8 ratel1, ratef1, ratel2, ratef2
      real*8 dnrm2, cnrml, cnrmf, ttgemml, ttgemmf
      external second_linux, dnrm2
c
c.....initialize arrays to zero
      do 10 i=1,nmax*kmax
         a(i) = one/float(i)
         c(i) = one/float(i)
         c1(i) = one/float(i)
 10   continue
      do 15 i=1,kmax*kmax
         b(i) = one/float(i)
         b1(i) = one/float(i)
 15   continue

      trans = 'Transpose'
      notrans = 'No Transpose'
c     
c.....Compute C = A*B
      beta = zero
      alpha = one
      n0 = 10000
      ntry = 2
      do 20 n1 = 1, 5
     
         if (n1.gt.1) then
            n0 = n0*2
         endif
         if (n0.gt.nmax) goto 20
         k0 = 1

         if (n1.eq.1) then
c     
c.....print header
            write(*,100)
     &           ' N ',' K ','BDOT ','BDOTF ','BAXPY ','BAXPYF '
            write(*,100)
     &           '---','---','-----','------','------','-------'
         endif

         do 25 k1 = 1,5
            if (k1.gt.1) then 
               k0 = k0*2
            endif
            if (k0.gt.kmax) goto 25
            if (k0.lt.4.and.n0.lt.16000) ntry = 4
            if (k0.lt.4.and.n0.lt.4000) ntry = 12
           
C     First compute (k by k) = (k by n) * (n by k) (Block dot product)
c     
c.....define number of ops to compute C = A*B
            fnops = 2.0*float(k0)*float(k0)*float(n0) 
     &              - float(k0)*float(k0)
            
c     
c.......time gemm
c     GEMM computes B = A(trans)*C (library routine)
         call second_linux(t1)
         do itry = 1,ntry
            call dgemm
     &           (trans,notrans,k0,k0,n0,alpha,a,n0,c,n0,beta,b,k0)
         end do
         call second_linux(t2)
         ttgemml = t2 - t1
         cnrml = dnrm2(m*n, b, 1)
c     
c.......time GEMMF
c     GEMMF computes C = beta*c + alpha*A*B (Fortran routine)
         call second_linux(t1)
         do itry = 1,ntry
            call dgemmf
     &           (trans,notrans,k0,k0,n0,alpha,a,n0,c,n0,beta,b1,k0)
         end do
         call second_linux(t2)
         ttgemmf = t2 - t1
         cnrmf = dnrm2(m*n, b1, 1)
c     
c.......Check answers
         if (abs(cnrml-cnrmf).gt.1.0D-08) then
            write(*,*)'WARN: Norms of library and Fortran differ'
            write(*,*)'WARN: Library = ',cnrml
            write(*,*)'WARN: Library = ',cnrmf
         endif
c     
c.....time for GEMM
         ratel1 = (float(ntry)*fnops)/(ttgemml*onemill)
c     
c.....time for GEMMF
         ratef1 = (float(ntry)*fnops)/(ttgemmf*onemill)
            
C     Second compute (n by k) = (n by k) * (k by k) (Block dot product)
c     
c.....define number of ops to compute C = A*B
            fnops = 2.0*float(n0)*float(k0)*float(k0) 
     &              - float(n0)*float(k0)
            
c     
c.......time gemm
c     GEMM computes C = A*B (library routine)
         call second_linux(t1)
         do itry = 1,ntry
            call dgemm
     &           (notrans,notrans,n0,k0,k0,alpha,a,n0,b,k0,beta,c,n0)
         end do
         call second_linux(t2)
         ttgemml = t2 - t1
         cnrml = dnrm2(m*n, b, 1)
c     
c.......time GEMMF
c     GEMMF computes C = beta*c + alpha*A*B (Fortran routine)
         call second_linux(t1)
         do itry = 1,ntry
            call dgemmf
     &           (notrans,notrans,n0,k0,k0,alpha,a,n0,b,k0,beta,c,n0)
         end do
         call second_linux(t2)
         ttgemmf = t2 - t1
         cnrmf = dnrm2(m*n, b1, 1)
c     
c.......Check answers
         if (abs(cnrml-cnrmf).gt.1.0D-08) then
            write(*,*)'WARN: Norms of library and Fortran differ'
            write(*,*)'WARN: Library = ',cnrml
            write(*,*)'WARN: Library = ',cnrmf
         endif
c     
c.....time for GEMM
         ratel2 = (float(ntry)*fnops)/(ttgemml*onemill)
c     
c.....time for GEMMF
         ratef2 = (float(ntry)*fnops)/(ttgemmf*onemill)
c     
c.....report results
         write(*,101)n0,k0,ratel1,ratef1,ratel2,ratef2
c     
c.....repeat
 25   continue
 20   continue
c     
 100  format(' ',a8,a8,a9,a9,1x,a7,1x,a7)
 101  format('  ',i8,i8,f8.2,1x,f8.2,1x,f8.2,f8.2)
 9999 continue
      end
