
c  subroutine icostr
c
      subroutine icostr(N, Fu, u, mmax, cosmat)
      implicit none
      integer  k, i, mmax
      double precision Fu, u, cosmat
      integer N 
      dimension Fu(0:mmax), u(0:mmax), cosmat(0:mmax, 0:mmax)
      double precision sum, PI
      PI=3.14159265358979d0

      do i=0, N
        sum=0.d0
        do k=0, N
          sum=sum + Fu(k)*cosmat(i,k)
        enddo
        u(i)=sum
      enddo
      return
      end

