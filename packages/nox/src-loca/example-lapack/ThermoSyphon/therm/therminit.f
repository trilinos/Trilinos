      subroutine therminit(res,jac,mass,m,mm,phi,psi,w,
     &               matR,matB,matI,matB2,matB2R,part1,
     &     lu_in, prandtl,rayleigh,bval1,bval2,file1,file2)
      implicit none
      include 'limits.h'
      include 'indata.h'
      include 'array.h'
      integer lu_in
      double precision rayleighread, prandtlread
      double precision part1(0:m_max, 0:m_max)
      double precision t, junk
      integer i, j
      integer res,jac,mass
      integer mm, m1
      character*16 file1,file2
      character*1 junkc

c Begin executables
      res=1
      jac=2
      mass=3
      open(unit=lu_in,file=file1,status='old',form='formatted')
        read(lu_in,*) junk 
        read(lu_in,*) m
        read(lu_in,*) junk
        read(lu_in,*) prandtl
        read(lu_in,*) junk 
        read(lu_in,*) rayleigh
        read(lu_in,*) bval1
        read(lu_in,*) bval2
        read(lu_in,*) junk   
        read(lu_in,*) junk   
        read(lu_in,*) junk   
        read(lu_in,*) junk    
        read(lu_in,*) junk     
        read(lu_in,*) junkc     
        read(lu_in,*) junkc     
      close(lu_in)
c     call read_indata(lu_in, m_max, n_max,m,prandtl,rayleigh,bval1,
c    &    bval2,file1,ra_step,ra_final,pr_final,jobvL,jobvR,tol)
c Initialize fields
      mm=m/2-1
      open(unit=lu_in, file=file2, status = 'unknown')
      read(lu_in,*) t
      read(lu_in,*) m1
      read(lu_in,*) rayleighread
      read(lu_in,*) prandtlread
      do i=0, min(m,m1)
        read(lu_in,*) phi(i), psi(i), w(i)
      enddo
      close(lu_in)
      do i=min(m,m1)+1, m
        phi(i)=0.d0
        psi(i)=0.d0
        w(i)=0.d0
      enddo
c     do i=0, m 
c       u(i)=phi(i)
c	u(i+m+1)=psi(i)
c       u(i+2*m+2)=w(i)
c     enddo
      do i=m+1, 2*m
        phi(i)=0.d0
        psi(i)=0.d0
        w(i)=0.d0
      enddo
C ends reading in
      call compR(m+buff,matR, m_max)
      call compB(m+buff,matB, m_max)
      call compI(m+buff,matI, m_max)
      call ffmul(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB,matB,matB2)
       do i=0, 1
        do j=0, m+buff
	 matB2(i,j)=0.d0
	enddo
       enddo
      call ffmul(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB2,matR,matB2R)
        do j=0, m+buff
          matR(0,j)=0.d0
          matR(1,j)=0.d0
          matB(0,j)=0.d0
          matB(1,j)=0.d0
        enddo
       do i=0, m
         do j=0, m
            part1(i,j)=matR(i,j)-matB(i,j)
         enddo
       enddo
      end




