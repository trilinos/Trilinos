C-----------------------------------------------------------------------
           SUBROUTINE FFMUL(MLIM,ML,KL,NL,A,B,C)
C-----------------------------------------------------------------------
C --- Multiply two matrices using a GAXPY algorithm (Gollub/VanLoan p.12)
C --- Assumes common leading column dimension (MLIM)
            IMPLICIT NONE
            INTEGER MLIM, ML, NL, KL
            DOUBLE PRECISION A(MLIM,KL)
            DOUBLE PRECISION B(MLIM,NL)
            DOUBLE PRECISION C(MLIM,NL)
            INTEGER I, J, K
C
            DO J = 1, NL
               DO I = 1, ML
                  C(I,J) = 0.D0
               END DO
            END DO
            DO J = 1, NL
               DO K = 1, KL
                  DO I = 1, ML
                     C(I,J) = C(I,J)+A(I,K)*B(K,J)
                  END DO
               END DO   
            END DO
           RETURN
           END
C-----------------------------------------------------------------------
      subroutine copym(m1,m2,m,n,indata,outdata)
C-----------------------------------------------------------------------
C Copy array into another; assumes different leading column dimensions
        integer   m1,m2,m,n
        double precision  indata(m1,*)
        DOUBLE PRECISION outdata(m2,*)
        integer     i, j
      do i=1, m
        do j = 1,n
          outdata(i,j) = indata(i,j)
        end do
      end do
      end      
C-----------------------------------------------------------------------
      subroutine mug(data,mmax,imin,imax,jmin,jmax)
C-----------------------------------------------------------------------
C-- breaks with the convention of this general-purpose package and 
C   dimensions arays from 0.
        implicit none
        integer i,j,imax,jmax,mmax,imin, jmin
        double precision data(0:mmax,0:*)
        do i = imin, imax
           write(*,*) i,(sngl(data(i,j)), j=jmin,jmax)
        end do
c1    format(I4,1x,F8.4)
      end
C-----------------------------------------------------------------------
      subroutine remmod(m_max,n_max,data,m_act,n_act)
C-----------------------------------------------------------------------
      implicit none
      integer m_max, n_max
      integer m_act, n_act
      double precision data(0:m_max,0:*)
      integer i, j
      do j = 0, n_act
         do i = m_act, m_max
            data(i,j) = 0.0d0
         end do
      end do
      do j = n_act, n_max+1
         do i = 0, m_max
            data(i,j) = 0.0d0
         end do
      end do
      do i = 0, m_act
         data(i,2) = 0.0d0
      end do
      end
C-----------------------------------------------------------------------
           subroutine point_mult(A,B,C,m_max,m,n)
C-----------------------------------------------------------------------
        integer m_max,m,n
        double precision A(0:m_max,0:*)
        double precision B(0:m_max,0:*)
        double precision C(0:m_max,0:*)
        integer i, j
        do i = 0,m
          do j = 0, n
             C(i,j) = A(i,j) * B(i,j)
          end do
        end do
        return
        end
C-----------------------------------------------------------------------
           subroutine addm(a0,A,b0,B,C,m_max,m,n)
C-----------------------------------------------------------------------
        integer m_max,m,n
        double precision A(0:m_max,0:*)
        double precision B(0:m_max,0:*)
        double precision C(0:m_max,0:*)
        double precision a0, b0
        integer i, j
        do j = 0, n
          do i = 0,m
             C(i,j) = a0*A(i,j) + b0*B(i,j)
          end do
        end do
        return
        end
C-----------------------------------------------------------------------
	subroutine error(mmax,m,n,xcomp,x,p,norm)
C-----------------------------------------------------------------------
C--compute a specified norm of the difference of two fields
C-----------------------------------------------------------------------
        implicit none
        integer mmax,m,n, p, i, j
	double precision xcomp(0:mmax,*), x(0:mmax,*), norm(*)
        do j = 0, n
          norm(j) = 0.d0
          if (p .eq. -1) then
            do i = 0,m
              norm(j)   = max(abs(xcomp(i,j)-x(i,j)),norm(j))
            end do
          else
            do i = 0,m
              norm(j) = norm(j) + (dabs(xcomp(i,j)-x(i,j)))**p
            end do
            norm(j) = norm(j)**(1.d0/dble(p))
          end if
          write(*,*) 'mode ',j,' error in L-',p,' norm: ',norm(j)
        end do
        return
        end
C-----------------------------------------------------------------------
