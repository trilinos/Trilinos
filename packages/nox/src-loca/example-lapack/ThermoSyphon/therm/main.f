      program newton
      implicit none
      integer*4 lu_in
      parameter(lu_in=17)
      integer*4 lu_old
      parameter(lu_old=19)
      integer*4 lu_print
      parameter(lu_print=15)
      integer*4 lu_time
      parameter(lu_time=21)
      integer*4 lu_name
      parameter (lu_name=11)
      integer*4 lu_e
      parameter (lu_e=3)
      include 'limits.h'
      include 'indata.h'
      include 'array.h'
      integer largestindex, li1, li2
      double precision eig,eig2,eig3,junk
      double precision rayleighread, prandtlread
      double precision work3(m_max)
      double precision evecphi(0:m_max)
      double precision evecpsi(0:m_max)
      double precision evecw(0:m_max)
      double precision evecr(3*m_max+3)
      double precision eveci(3*m_max+3)
      double precision cosmat(0:m_max, 0:m_max)
      double precision packOp(0:3*m_max+2, 0:3*m_max+2)
      double precision part1(0:m_max, 0:m_max)
      double precision convw(0:m_max, 0:m_max)
      double precision convpsi(0:m_max, 0:m_max)
      double precision convphi(0:m_max, 0:m_max)
      double precision block12(0:m_max, 0:m_max)
      double precision block13(0:m_max, 0:m_max)
      double precision block23(0:m_max, 0:m_max)
      double precision block31(0:m_max, 0:m_max)
      double precision psiw(0:m_max)
      double precision phiw(0:m_max)
      double precision lapphi(0:m_max)
      double precision lappsi(0:m_max)
      double precision lapw(0:m_max)
      double precision rhs1(0:m_max)
      double precision rhs2(0:m_max)
      double precision rhs3(0:m_max)
      double precision x(0:m_max)
      double precision rhs(0:3*m_max+2)
      double precision packrhs(0:3*m_max+2)
      double precision spacephi(0:m_max)
      double precision spacepsi(0:m_max)
      double precision spacew(0:m_max)
      double precision u(0:3*m_max+2)
      double precision packu(0:3*m_max+2)
      double precision sol(0:3*m_max+2)
      double precision  work(3*m_max+2), rcond
      double precision actual_error
      double precision t
      integer iter
      integer ipvt(3*m_max+2+1), job
      integer i, j
      integer res,jac,mass
      integer mm, m1
      double precision matM(0:3*m_max+2, 0:3*m_max+2)
      double precision matEig(0:3*m_max+2, 0:3*m_max+2)
      double precision  work2(8*(3*m_max+3))
      double precision vL(3*m_max+3,3*m_max+3),vR(3*m_max+3,3*m_max+3)
      integer ldvL, ldvR
      integer lwork, info
      double precision alphai(3*m_max+2+1), alphar(3*m_max+2+1)
      double precision beta(3*m_max+2+1)
      character*16 file1,file2

      

c Begin executables
      lwork=8*(3*m_max+3)
      ldvl=3*m_max+3
      ldvr=3*m_max+3
      file1='in'
      file2='oldfile'
      call therminit(res,jac,mass,m,mm,phi,psi,w,
     &               matR,matB,matI,matB2,matB2R,part1,
     &     lu_in, prandtl,rayleigh,bval1,bval2,file1,file2)
        open(unit=lu_in,file=file1,status='old',form='formatted')
          read(lu_in,*) debug
          read(lu_in,*) junk
          read(lu_in,*) n
          read(lu_in,*) junk
          read(lu_in,*) sigma
          read(lu_in,*) junk    
          read(lu_in,*) junk    
          read(lu_in,*) junk    
          read(lu_in,*) bc_case
          read(lu_in,*) tol
          read(lu_in,*) ra_step
          read(lu_in,*) ra_final
          read(lu_in,*) pr_final
          read(lu_in,*) jobvL
          read(lu_in,*) jobvR
        close(lu_in)
      do i=0, m
        u(i)=phi(i)
        u(i+m+1)=psi(i)
        u(i+2*m+2)=w(i)
      enddo
      call initcos(m, cosmat, m_max)
      call icostr(m, w, spacew, m_max, cosmat)
      call icostr(m, phi, spacephi, m_max, cosmat)
      call icostr(m, psi, spacepsi, m_max, cosmat)
c ends making building blocks
c outer loop
5     if (rayleigh .gt. ra_final) goto 30
       print*, 'ra', rayleigh, 'ra_step', ra_step, 'rafinal', ra_final
      actual_error=1.0
      iter=0
c now enter the loop
10    if ((actual_error .lt. tol) .or. (iter .gt.10) ) goto 20
      do i=0,mm
         packu(i)= phi(2*i)
         packu(i+mm+1)=psi(2*i) 
         packu(i+2*mm+2)=w(2*i) 
      enddo
      call thermfill(res,m_max,mm,packu,w,phi,psi,convw,convphi,
     &        convpsi,matB2R,rayleigh,prandtl,block12,block13,
     &        block23,block31,part1,lapphi,lappsi,lapw,phiw,psiw,
     &        matB,matR,matI,
     &        rhs1,rhs2,rhs3,bval1,bval2,m,buff,
     &        packOp,packrhs,matM)
      do i=0, 3*mm+2
         print*, packu(i), packrhs(i)
      enddo
      call thermfill(jac,m_max,mm,packu,w,phi,psi,convw,convphi,
     &        convpsi,matB2R,rayleigh,prandtl,block12,block13,
     &        block23,block31,part1,lapphi,lappsi,lapw,phiw,psiw,
     &        matB,matR,matI,
     &        rhs1,rhs2,rhs3,bval1,bval2,m,buff,
     &        packOp,packrhs,matM)
c      do i=0,mm
c         packrhs(i)=resphi(2*i)
c         packrhs(i+mm+1)=respsi(2*i)
c         packrhs(i+2*mm+2)=resw(2*i)
c      enddo
      call thermfill(mass,m_max,mm,packu,w,phi,psi,convw,convphi,
     &        convpsi,matB2R,rayleigh,prandtl,block12,block13,
     &        block23,block31,part1,lapphi,lappsi,lapw,phiw,psiw,
     &        matB,matR,matI,
     &        rhs1,rhs2,rhs3,bval1,bval2,m,buff,
     &        packOp,packrhs,matM)
c     do i=0, mm
c       do j=0,mm
c        matM(i,j)=unpmatM(2*i+1,2*j)
c        matM(i+mm+1, j+mm+1) = unpmatM(2*(i+mm+1)+1,2*(j+mm+1))
c        matM(i+2*mm+2,j+2*mm+2) = unpmatM(2*(i+2*mm+2)+1,2*(j+2*mm+2))
c       enddo
c      enddo
c boundary condition
c       do j=0, mm
c         matM(0,j)=0.d0
c         matM(mm+1, j+mm+1)=0.d0
c         matM(2*mm+2, j+2*mm+2)=0.d0
c       enddo


       do i=0, 3*mm+2
         do j=0, 3*mm+2
            matEig(i,j)=packOp(i,j)
         enddo
       enddo
c factor
      job=0
c      call dgeco(packOp,3*m_max+2+1,3*mm+2+1,ipvt,rcond,work)
c solve
c      call dgesl(packOp,3*m_max+2+1,3*mm+2+1,ipvt,packrhs ,job)
      call dgesv(3*mm+3,1,packOp,3*m_max+3,ipvt,packrhs,3*m_max+3,job)
      actual_error=abs(packrhs(0))
      do i=1,3*mm+2
       if (abs(packrhs(i)) .gt. actual_error) then
         actual_error = abs(packrhs(i))
       endif
      enddo
      print*, 'act error', actual_error
c unpack packrhs
      do i=0, 3*m+2
        rhs(i)=0.d0
      enddo
      do i=0, mm
        rhs(2*i)=packrhs(i)
        rhs(2*i+1+m)=packrhs(i+mm+1)
        rhs(2*i+2+2*m)=packrhs(i+2*mm+2)
      enddo
      do i=0, 3*m+2
        u(i)=u(i)+rhs(i)
      enddo
      
c fills the ends with zeros for convol. accuracy
      do i=0,2*m
       phi(i)=0.d0
       psi(i)=0.d0
       w(i)=0.d0
      enddo
      do i=0, m
        phi(i)=u(i)
        psi(i)=u(i+m+1)
        w(i)=u(i+2*m+2)
c       print*, i, phi(i), psi(i),  iter
      enddo
      iter=iter+1
      goto 10
20    continue
c print out solution in modespace (for reading in in timeint code)
      open(unit=lu_time, file='oldfile', status = 'unknown')
      write(lu_time, *) 0
      write(lu_time,*) m
      write(lu_time,*) rayleigh
      write(lu_time,*) prandtl 
      do i=0,m
        write(lu_time,*) phi(i)
        write(lu_time,*) psi(i)
        write(lu_time,*) w(i)
      enddo
      open(unit=8,file='chebcoef', status='unknown')
      do i=0,m
        write(8,*) phi(i)
        write(8,*) psi(i)
      enddo
      close(17)
      call icostr(m, w, spacew, m_max, cosmat)
      call icostr(m, phi, spacephi, m_max, cosmat)
      call icostr(m, psi, spacepsi, m_max, cosmat)
      open(unit=3, file='vel.m')
      open(unit=4, file='temp.m')
c     do i=0, m
c      write(3,*) 'V(', i+1, ',:)=[', x(i), spacew(i), '];'
c     write(4,*) 'T(', i+1, ',:)=[', spacephi(i), spacepsi(i),'];'
c     enddo
c some file printing stuff
c     open(unit=lu_name, file='namefile', status='unknown')
c     if (Rayleigh .gt. 100) then
c     write(unit=lu_name, fmt=41) prandtl, Rayleigh
c     elseif( (Rayleigh .ge.10).and.(Rayleigh .lt.100) ) then
c     write(unit=lu_name, fmt=42) prandtl, Rayleigh
c     endif
41    format('vprofpr', f5.3, 'ra', f6.2)
42    format('vprofpr', f5.3, 'ra', f5.2)
c     close(unit=lu_name)
c     open(unit=lu_name, file='namefile', status='old')
c     read(lu_name,*) graphfile
c     close(unit=lu_name)
c     do i=0,m
c      print*, rayleigh, spacepsi(i) 
c     enddo
c      write(lu_name,*)rayleigh, spacepsi(lind-1) 
c      print*, rayleigh, lind
c     close(unit=lu_print)
      rayleigh=rayleigh+1.d0*ra_step
c do eigenvalue problem
      open(unit=lu_print, file='eigmats', status = 'unknown')
      write(lu_print,*) mm+1
      do i=0,3*mm+2
        do j=0, 3*mm+2
          write(lu_print,*)  matEig(i,j)
	enddo
      enddo
      do i=0,3*mm+2
        do j=0, 3*mm+2
          write(lu_print,*)  matM(i,j)
	enddo
      enddo
      close(unit=lu_print)
      info=0
      call dggev(jobvL,jobvR,3*mm+2+1,matEig,3*m_max+2+1,matM,
     &            3*m_max+2+1,alphar, alphai, beta, vL, ldvL, 
     &            vR, ldvR, work2, lwork, info)
       print*, 'info', info
       do i=1, 3*mm+3
         if(beta(i) .ne. 0.d0) then
	  work3(i)=alphar(i)/beta(i)
	 else
c set it very small so won't interfere with sort
           work3(i)=-100000
	 endif
       enddo
       call findlargest(work3, 3*mm+3, largestindex)
       print*, 'eigs'
      do i=1, 3*mm+3
         print*,i, alphar(i)/beta(i), alphai(i)/beta(i)
c    &     , beta(i)
      enddo
	  eig=alphar(largestindex)/beta(largestindex)
	  eig2=alphar(largestindex-1)/beta(largestindex-1)
	  eig3=alphar(largestindex+1)/beta(largestindex+1)
       if (abs(eig-eig2) .lt. tol) then
          li1=largestindex-1
	  li2=largestindex
       elseif (abs(eig-eig3) .lt. tol) then
           li1=largestindex
	   li2=largestindex+1
       else 
          print*, 'wrong eval index'
       endif
         do i=1,3*mm+3
	   evecr(i)=vR(i,li1)
	   eveci(i)=vR(i,li2)
	 enddo
c print out right eigenvector in modespace (for reading in timeint code)
c extract the phi,psi, w parts, imaginary part
      do i=0,m
        evecphi(i)=0.d0
        evecpsi(i)=0.d0
        evecw(i)=0.d0
      enddo
      do i=0, mm
        evecphi(2*i)=vR(i+1,li2)
        evecpsi(2*i)=vR(i+1+mm+1,li2)
        evecw(2*i)=vR(i+1+2*mm+2,li2)
      enddo
      open(unit=lu_e, file='evecs', status = 'unknown')
      write(lu_e, *) 0
      write(lu_e,*) m
      write(lu_e,*) rayleigh-ra_step
      write(lu_e,*) prandtl 
      write(lu_e,*) alphai(li1)/beta(li1) 
      do i=0,m
        write(lu_e,*) evecphi(i)
        write(lu_e,*) evecpsi(i)
        write(lu_e,*) evecw(i)
      enddo
      close(lu_e)
      open(unit=lu_e, file='omega', status = 'unknown')
      write(lu_e,*) alphai(li1)/beta(li1) 
      close(lu_e)
      goto 5
30    continue

      end




