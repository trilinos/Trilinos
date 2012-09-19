      subroutine asol_ug3d_so_ratedef( nelem,nint,dt,
     *     def_grad_inv, cg_left_inv, eval, evec, rate_def )
c
c***********************************************************************
c***********************************************************************
c
c description:
c     This routine computes the rate of deformation of the material using
c     a "strong objectivity" approach.  More description to follow.
c
c formal parameters - input:
c     nelem              int   number of elements in the workset
c     nint               int   number of integration stations for element
c     dt                 Real  time increment
c     def_grad_inv       Real  velocity gradient
c
c formal parameters - workspace for computations:
c     cg_left_inv        Real  inverse of left cauchy green tensor (B^-1)
c     eval               Real  eigenvalues of B^-1
c     evec               Real  eigenvectors of B^-1
c     stretch_left_inv   Real  inverse of left stretch (V^-1)
c
c formal parameters - output:
c     def_grad_inv       Real  inverse of incremental deformation gradient
c     rate_def           Real  rate of deformation (D)
c
c***********************************************************************
c
#include <precision.par>
#include <numbers.par>
#include <framewk/Fmwk_type_sizes.par>
#include <math/Fmth_math.par>
c
c         mappings from sierra structure to local structure
c
      parameter( iv_x = k_v3dx,
     2           iv_y = k_v3dy,
     2           iv_z = k_v3dz )
      parameter( is_xx = k_s33xx,
     2           is_yy = k_s33yy,
     2           is_zz = k_s33zz,
     2           is_xy = k_s33xy,
     2           is_yz = k_s33yz,
     2           is_zx = k_s33zx )
      parameter( if_xx = k_f36xx,
     *           if_yy = k_f36yy,
     *           if_zz = k_f36zz,
     *           if_xy = k_f36xy,
     *           if_yz = k_f36yz,
     *           if_zx = k_f36zx,
     *           if_yx = k_f36yx,
     *           if_zy = k_f36zy,
     *           if_xz = k_f36xz )
c
c         input parameters
c
      dimension def_grad_inv(n_f36,nelem,nint),
     *     cg_left_inv(n_s33,nelem,nint),
     *     eval(n_v3d,nelem,nint),evec(n_f36,nelem,nint),
     *     rate_def(n_s33,nelem,nint)
c
c        Convert velocity gradient into incremental rate of deformation
c           ( F^-1 = I - du/dx )
c
      do j = 1, nint
         do k = 1, nelem
c
            def_grad_inv(if_xx,k,j)=-def_grad_inv(if_xx,k,j)*dt
     &           + ONE
            def_grad_inv(if_yy,k,j)=-def_grad_inv(if_yy,k,j)*dt
     &           + ONE
            def_grad_inv(if_zz,k,j)=-def_grad_inv(if_zz,k,j)*dt
     &           + ONE
            def_grad_inv(if_xy,k,j)=-def_grad_inv(if_xy,k,j)*dt
            def_grad_inv(if_yz,k,j)=-def_grad_inv(if_yz,k,j)*dt
            def_grad_inv(if_zx,k,j)=-def_grad_inv(if_zx,k,j)*dt
            def_grad_inv(if_yx,k,j)=-def_grad_inv(if_yx,k,j)*dt
            def_grad_inv(if_zy,k,j)=-def_grad_inv(if_zy,k,j)*dt
            def_grad_inv(if_xz,k,j)=-def_grad_inv(if_xz,k,j)*dt
c
         enddo
      enddo
c
c        If dt is zero, then there is no rate of deformation for this
c        loading increment. Set rate of deformation tensor to be zero.
c
      if ( dt .eq. ZERO) then
c
         do j = 1, nint
            do k = 1, nelem
c
               rate_def(is_xx,k,j) = ZERO
               rate_def(is_yy,k,j) = ZERO
               rate_def(is_zz,k,j) = ZERO
               rate_def(is_xy,k,j) = ZERO
               rate_def(is_yz,k,j) = ZERO
               rate_def(is_zx,k,j) = ZERO
c
            enddo
         enddo
         return
c
      endif
c
c        Compute inverse of left cauchy green tensor (B^-1)
c
      call fmath_right_tensor_square36( nelem,nint,def_grad_inv,
     &     cg_left_inv)
c
c        Get eigenvalues of B^-1
c
      call fmth_eigen( nelem,nint,cg_left_inv,eval,evec )
c
      do j = 1, nint
         do k = 1, nelem
c
            xlog = log(eval(iv_x,k,j))
            ylog = log(eval(iv_y,k,j))
            zlog = log(eval(iv_z,k,j))
c
            xinv = sqrt(eval(iv_x,k,j))
            yinv = sqrt(eval(iv_y,k,j))
            zinv = sqrt(eval(iv_z,k,j))
c
            rate_def(is_xx,k,j) = - (half/dt)*(
     &       xlog*evec(if_xx,k,j)*evec(if_xx,k,j)
     &      +ylog*evec(if_xy,k,j)*evec(if_xy,k,j)
     &      +zlog*evec(if_xz,k,j)*evec(if_xz,k,j))
c
            rate_def(is_yy,k,j) = - (half/dt)*(
     &       xlog*evec(if_yx,k,j)*evec(if_yx,k,j)
     &      +ylog*evec(if_yy,k,j)*evec(if_yy,k,j)
     &      +zlog*evec(if_yz,k,j)*evec(if_yz,k,j))
c
            rate_def(is_zz,k,j) = - (half/dt)*(
     &       xlog*evec(if_zx,k,j)*evec(if_zx,k,j)
     &      +ylog*evec(if_zy,k,j)*evec(if_zy,k,j)
     &      +zlog*evec(if_zz,k,j)*evec(if_zz,k,j))
c
            rate_def(is_xy,k,j) = - (half/dt)*(
     &       xlog*evec(if_xx,k,j)*evec(if_yx,k,j)
     &      +ylog*evec(if_xy,k,j)*evec(if_yy,k,j)
     &      +zlog*evec(if_xz,k,j)*evec(if_yz,k,j))
c
            rate_def(is_yz,k,j) = - (half/dt)*(
     &       xlog*evec(if_yx,k,j)*evec(if_zx,k,j)
     &      +ylog*evec(if_yy,k,j)*evec(if_zy,k,j)
     &      +zlog*evec(if_yz,k,j)*evec(if_zz,k,j))
c
            rate_def(is_zx,k,j) = - (half/dt)*(
     &       xlog*evec(if_zx,k,j)*evec(if_xx,k,j)
     &      +ylog*evec(if_zy,k,j)*evec(if_xy,k,j)
     &      +zlog*evec(if_zz,k,j)*evec(if_xz,k,j))
c
         enddo
      enddo
c
      return
      end
c
c     Copyright 2000 Sandia Corporation, Albuquerque, NM.
c
