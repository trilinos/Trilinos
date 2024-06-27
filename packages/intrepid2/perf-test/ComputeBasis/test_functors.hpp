// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_functors.cpp
    \brief  Performance test comparing dynrankview overhead
    \author Created by Kyungjoo Kim.
*/

namespace Intrepid2 {
  
  namespace Test {

    template<typename VectorType,
             typename ValueType,
             typename DeviceSpaceType>
    struct F_hgrad_eval {
      typedef Kokkos::DynRankView<VectorType,DeviceSpaceType> VectorViewType;
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> ViewType;
      typedef typename ViewType::const_type constViewType;

      typedef VectorType vector_type;
      typedef ValueType value_type;

      /**/  VectorViewType _weighted_basis_values;
      /**/  VectorViewType _weighted_basis_grads;
      const       ViewType _grads;
      const VectorViewType _workset;
      const       ViewType _weights;
      const       ViewType _basis_values;
      const       ViewType _basis_grads;

      KOKKOS_INLINE_FUNCTION
      F_hgrad_eval(/**/  VectorViewType weighted_basis_values_,
                   /**/  VectorViewType weighted_basis_grads_,
                   const       ViewType grads_, // gradient for jacobian
                   const VectorViewType workset_, // workset input
                   const       ViewType weights_, // weights
                   const       ViewType basis_values_, // reference values
                   const       ViewType basis_grads_) // reference grads
        : _weighted_basis_values(weighted_basis_values_),
          _weighted_basis_grads(weighted_basis_grads_),
          _grads(grads_),
          _workset(workset_),
          _weights(weights_),
          _basis_values(basis_values_),
          _basis_grads(basis_grads_) {}

      KOKKOS_INLINE_FUNCTION
      void apply(const ordinal_type cl,
                 const ordinal_type pt) const {
        vector_type buf[9 + 9];

        const auto grad = Kokkos::subview(_grads,       Kokkos::ALL(), pt, Kokkos::ALL());
        const auto dofs = Kokkos::subview(_workset, cl, Kokkos::ALL(),     Kokkos::ALL());

        const ordinal_type card = dofs.extent(0);
        const ordinal_type dim = dofs.extent(1);
        
        // temporary values
        Kokkos::View<vector_type**, Kokkos::AnonymousSpace> 
          jac    (&buf[0], dim, dim), 
          jac_inv(&buf[9], dim, dim); 

        // setJacobian  F_setJacobian::apply(jac, dofs, grads);
        {
          for (ordinal_type i=0;i<dim;++i)
            for (ordinal_type j=0;j<dim;++j) {
              jac(i, j) = 0;
              for (ordinal_type bf=0;bf<card;++bf)
                jac(i, j) += dofs(bf, i)*grad(bf, j);
            }
        }
        
        // setJacobianDet   F_setJacobianDet::apply(det, jac);
        // setJacobianInv  F_setJacobianInv::apply(jac_inv, jac);
        vector_type det = 0;
        if (dim == 2) {
          det = ( jac(0,0) * jac(1,1) -
                  jac(0,1) * jac(1,0) );
        } else {
          det = ( jac(0,0) * jac(1,1) * jac(2,2) +
                  jac(1,0) * jac(2,1) * jac(0,2) +
                  jac(2,0) * jac(0,1) * jac(1,2) -
                  jac(2,0) * jac(1,1) * jac(0,2) -
                  jac(0,0) * jac(2,1) * jac(1,2) -
                  jac(1,0) * jac(0,1) * jac(2,2) );
        }

        if (dim == 2) {
          const vector_type val = det;
          jac_inv(0,0) =   jac(1,1)/val;
          jac_inv(1,1) =   jac(0,0)/val;

          jac_inv(1,0) = -1.0* jac(1,0)/val;
          jac_inv(0,1) = -1.0* jac(0,1)/val;
        } else {
          const vector_type val = det;
          vector_type val0, val1, val2;

          val0 =   jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2);
          val1 =   jac(2,0)*jac(1,2) - jac(1,0)*jac(2,2);
          val2 =   jac(1,0)*jac(2,1) - jac(2,0)*jac(1,1);

          jac_inv(0,0) = val0/val;
          jac_inv(1,0) = val1/val;
          jac_inv(2,0) = val2/val;

          val0 =   jac(2,1)*jac(0,2) - jac(0,1)*jac(2,2);
          val1 =   jac(0,0)*jac(2,2) - jac(2,0)*jac(0,2);
          val2 =   jac(2,0)*jac(0,1) - jac(0,0)*jac(2,1);

          jac_inv(0,1) = val0/val;
          jac_inv(1,1) = val1/val;
          jac_inv(2,1) = val2/val;

          val0 =   jac(0,1)*jac(1,2) - jac(1,1)*jac(0,2);
          val1 =   jac(1,0)*jac(0,2) - jac(0,0)*jac(1,2);
          val2 =   jac(0,0)*jac(1,1) - jac(1,0)*jac(0,1);

          jac_inv(0,2) = val0/val;
          jac_inv(1,2) = val1/val;
          jac_inv(2,2) = val2/val;
        }
        
        // computeCellMeasure
        const vector_type cell_measure = det*_weights(pt);
        
        // multiplyMeasure
        for (ordinal_type bf=0;bf<card;++bf)        
          _weighted_basis_values(cl, bf, pt) = cell_measure*_basis_values(bf, pt);
        
        // HGRADtransformGRAD
        auto weighted_grad = Kokkos::subview(_weighted_basis_grads, cl, Kokkos::ALL(), pt, Kokkos::ALL());
        const auto basis_grad = Kokkos::subview(_basis_grads, Kokkos::ALL(), pt, Kokkos::ALL());
        
        {
          const ordinal_type card = basis_grad.extent(0);
          const ordinal_type dim = basis_grad.extent(1);
          for (ordinal_type bf=0;bf<card;++bf)
            for (ordinal_type i=0;i<dim;++i) {
              weighted_grad(bf, i) = 0;
              for (ordinal_type j=0;j<dim;++j) 
                weighted_grad(bf, i) += jac_inv(i,j)*basis_grad(bf, j);
              weighted_grad(bf, i) *= cell_measure;
            }
        }
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl, 
                      const ordinal_type pt) const {
        apply(cl, pt);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type npts = _basis_values.extent(1);
        for (ordinal_type pt=0;pt<npts;++pt) 
          apply(cl, pt);
      }
    };
  } // end of namespace TEST
} // end of namespace Intrepid2
