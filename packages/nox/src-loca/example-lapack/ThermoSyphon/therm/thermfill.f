      subroutine thermfill(flag,m_max,mm,packu,w,phi,psi,convw,convphi,
     &        convpsi,matB2R,lambda,prandtl,block12,block13,
     &        block23,block31,part1,lapphi,lappsi,lapw,phiw,psiw,
     &        matB,matR,matI,
     &        rhs1,rhs2,rhs3,bval1,bval2,m,buff,
     &         packOp,packrhs,matM) 
      implicit none
      integer m_max
      include 'array.h'
      integer buff
      integer flag
      double precision lambda,prandtl,bval1,bval2
      integer m
      double precision packOp(0:3*m_max+2, 0:3*m_max+2)
      double precision packu(0:3*m_max+2)
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
      double precision packrhs(0:3*m_max+2)
      integer i, j
      integer mm
      double precision matM(0:3*m_max+2, 0:3*m_max+2)

c
c unpack
c
      do i=0, mm
        phi(2*i)=packu(i)
        psi(2*i)=packu(i+mm+1)
        w(2*i)=packu(i+2*mm+2)
      enddo
c create convolution
c convw holds convolution with w
      call convol(m_max,m,w,convw)
c convpsi holds convolution with psi
      call convol(m_max,m,psi,convpsi)
c convphi holds convolution with phi
      call convol(m_max,m,phi,convphi)
c create the blocks of the jacobian
c
      call FFMUL(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB2R,
     &       convw,block12) 
      call FFMUL(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB2R,
     &       convpsi,block13) 
      call FFMUL(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB2R,
     &       convphi,block23) 
      call FFMUL(m_max+1,m+1+buff,m+1+buff,m+1+buff,matB2R,
     &       matI,block31) 
      if (flag.eq.2) then
      do i=0,3*mm+2
       do j=0,3*mm+2
         packOp(i,j)=0.d0
       enddo
      enddo
       do i=0, mm
        do j=0, mm
          packOp(i,j)=matR(2*i+1,2*j)-matB(2*i+1,2*j)
          packOp(i, j+mm+1)=lambda*block12(2*i+1,2*j)
          packOp(i, j+2*mm+2)=lambda*block13(2*i+1,2*j)
          packOp(i+mm+1, j)=-lambda*block12(2*i+1,2*j)
          packOp(i+mm+1, j+mm+1)=matR(2*i+1,2*j)-matB(2*i+1,2*j)
          packOp(i+mm+1, j+2*mm+2)=-lambda*block23(2*i+1,2*j)
          packOp(i+2*mm+2, j)=prandtl*block31(2*i+1,2*j)
          packOp(i+2*mm+2,j+mm+1)=0.d0
          packOp(i+2*mm+2,j+2*mm+2)=matR(2*i+1,2*j)-matB(2*i+1,2*j)
        enddo
       enddo
c boundary condition
       do j=0, mm
         packOp(0,j)=1.d0
         packOp(mm+1, j+mm+1)=1.d0
         packOp(2*mm+2, j+2*mm+2)=1.d0
       enddo
       elseif (flag.eq.1) then
c create rhs      
c psiw holds psi*w
c phiw holds phi*w
       call ffmul(m_max+1, m+1, m+1,1, part1, phi, lapphi)
       call ffmul(m_max+1, m+1, m+1,1, part1, psi, lappsi)
       call ffmul(m_max+1, m+1, m+1,1, part1, w, lapw)
       call FFMUL(m_max+1,m+1,m+1,1,convpsi,w,psiw) 
       call FFMUL(m_max+1,m+1,m+1,1,convphi,w,phiw) 
       call FFMUL(m_max+1,m+1,m+1,1,matB2R, psiw, rhs1)
       call FFMUL(m_max+1,m+1,m+1,1,matB2R, phiw, rhs2)
       call FFMUL(m_max+1,m+1+buff,m+1+buff,1,matB2R, phi, rhs3)
c      do i=0,m
c         resphi(i)=-(lapphi(i+1)+lambda*rhs1(i+1))
c         respsi(i)=-(lappsi(i+1)-lambda*rhs2(i+1))
c         resw(i)=-(lapw(i+1)+prandtl*rhs3(i+1))
c      enddo
       do i=0, 3*m+2
        packrhs(i)=0.d0
       enddo
       do i=0,mm
         packrhs(i)=-(lapphi(2*i+1)+lambda*rhs1(2*i+1))
         packrhs(i+mm+1)=-(lappsi(2*i+1)-lambda*rhs2(2*i+1))
         packrhs(i+2*mm+2)=-(lapw(2*i+1)+prandtl*rhs3(2*i+1))
       enddo 
c boundary condition
c         resphi(0)=bval1
c         respsi(0)=bval2
c         resw(0)=bval1
       packrhs(0)=bval1
       packrhs(mm+1)=bval2
       packrhs(2*mm+2)=bval1
       elseif (flag.eq.3) then
c eigenvalue stuff
c create M matrix
      do i=0, 3*m_max+2
       do j=0, 3*m_max+2
         matM(i,j)=0.d0
       enddo
      enddo
c     do i=0,m
c       do j=0,m
c         matM(i,j)=prandtl*block31(i,j)
c         matM(i+m+1,j+m+1)=prandtl*block31(i,j)
c         matM(i+2*m+2,j+2*m+2)=block31(i,j)
c       enddo
c     enddo
c      do j=0, m
c        matM(0,j)=0.d0
c        matM(m+1, j+m+1)=0.d0
c        matM(2*m+2, j+2*m+2)=0.d0
c      enddo
c block31 holds B2R*I
       do i=0, mm
        do j=0,mm
         matM(i,j)=prandtl*block31(2*i+1,2*j)
         matM(i+mm+1, j+mm+1) = prandtl*block31(2*i+1,2*j)
         matM(i+2*mm+2,j+2*mm+2) = block31(2*i+1,2*j)
        enddo
       enddo
c boundary condition
        do j=0, mm
          matM(0,j)=0.d0
          matM(mm+1, j+mm+1)=0.d0
          matM(2*mm+2, j+2*mm+2)=0.d0
        enddo
      else
        print*, 'Bad value of fill flag'
      endif
      do i=0,mm
         packu(i)= phi(2*i)
         packu(i+mm+1)=psi(2*i)
         packu(i+2*mm+2)=w(2*i)
      enddo


      end





