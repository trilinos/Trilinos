c vec has to come in twice as long as the mat is wanted to avoid 
c "overlap" error
C-----------------------------------------------------------------------
        subroutine convol(mmax,m,vec,mat)
C-----------------------------------------------------------------------
        implicit none
        integer mmax, m, i, k
        double precision vec(0:2*mmax), mat(0:mmax,0:*)
        do i = 0,mmax
          do k = 0,mmax
            mat(i,k) = 0.d0
          end do
        end do
        do k = 0,m
          do i = k,m
            mat(i,i-k) = mat(i,i-k)+.5d0*vec(k)
            mat(i-k,i) = mat(i-k,i)+.5d0*vec(k)
          end do
          do i = 1, k
            mat(i,k-i) = mat(i,k-i)+.5d0*vec(k)
          end do
        end do
        do k = m+1, 2*m
          do i = k-m, m
            mat(i,k-i) = mat(i,k-i)+.5d0*vec(k)
          end do
        end do
        return
        end
C-----------------------------------------------------------------------

