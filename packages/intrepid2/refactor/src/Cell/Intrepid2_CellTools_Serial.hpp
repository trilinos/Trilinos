// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_CellTools_Serial.hpp
    \brief  Header file for the Impl::Intrepid2::CellTools::Serial class.
    \author Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOOLS_SERIAL_HPP__
#define __INTREPID2_CELLTOOLS_SERIAL_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Kernels.hpp"

namespace Intrepid2 {

  namespace Impl {
    
    class CellTools {
    public:
      struct Serial {

        // output: 
        //   jacobian (D,D) - jacobian matrix evaluated at a single point 
        // input: 
        //   grads    (N,D) - hgrad basis grad values evaluated at a single point (C1/C2 element only)
        //   nodes    (N,D) - cell element-to-node connectivity
        template<typename jacobianViewType,
                 typename basisGradViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        computeJacobian(const jacobianViewType  &jacobian, // D,D
                        const basisGradViewType &grads,    // N,D
                        const nodeViewType      &nodes) {  // N,D
          const auto N = nodes.dimension_0();
          const auto D = jacobian.dimension_0();
          
          INTREPID2_TEST_FOR_ABORT(jacobian.dimension_0() != jacobian.dimension_1(), "jacobian is not a square matrix.");
          INTREPID2_TEST_FOR_ABORT(N != grads.dimension_0(), "grad dimension_0 does not match to cardinality.");
          INTREPID2_TEST_FOR_ABORT(D != grads.dimension_1(), "grad dimension_1 does not match to space dim.");
          INTREPID2_TEST_FOR_ABORT(D != nodes.dimension_1(), "node dimension_1 does not match to space dim.");

          Kernels::Serial::gemm_trans_notrans(1.0, nodes, grads, 0.0, jacobian);
        }

        // output: 
        //   point (D)   - mapped physical point 
        // input: 
        //   vals  (N)   - hgrad basis values evaluated at a single point (C1/C2 element only)
        //   nodes (N,D) - cell element-to-node connectivity
        template<typename pointViewType,
                 typename basisValViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        mapToPhysicalFrame(const pointViewType    &point,    // D  
                           const basisValViewType &vals,     // N  
                           const nodeViewType     &nodes) {  // N,D 
          const auto N = vals.dimension_0();
          const auto D = point.dimension_0();

          INTREPID2_TEST_FOR_ABORT(N != nodes.dimension_0(), "nodes dimension_0 does not match to vals dimension_0.");
          INTREPID2_TEST_FOR_ABORT(D != nodes.dimension_1(), "node dimension_1 does not match to space dim.");

          Kernels::Serial::gemv_trans(1.0, nodes, vals, 0.0, point);
        }

        // template:
        //   implBasisType - impl basis function type e.g., Impl::Basis_HGRAD_QUAD_C1_FEM
        // output: 
        //   xref (D)    - point mapped to reference frame
        // input: 
        //   xphy  (D)   - point in physical frame
        //   nodes (N,D) - cell element-to-node connectivity
        template<typename implBasisType,
                 typename refPointViewType,
                 typename phyPointViewType,
                 typename nodeViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        mapToReferenceFrame(const refPointViewType &xref, // D
                            const phyPointViewType &xphy, // D
                            const nodeViewType &nodes) {  // N,D
          const ordinal_type N = nodes.dimension_0();
          const ordinal_type D = nodes.dimension_1();
          
          INTREPID2_TEST_FOR_ABORT(D != static_cast<ordinal_type>(xref.dimension_0()), "xref dimension_0 does not match to space dim.");
          INTREPID2_TEST_FOR_ABORT(D != static_cast<ordinal_type>(xphy.dimension_0()), "xphy dimension_0 does not match to space dim.");
          
          typedef typename refPointViewType::non_const_value_type value_type;
          
          // I want to use view instead of dynrankview
          // NMAX = 28, MAXDIM = 3
          value_type buf[27*3 + 27 + 9 + 9 + 3 + 3] = {}, *ptr = &buf[0];
          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> grads(ptr, N, D); ptr += N*D;
          
          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> vals(ptr, N); ptr += N;

          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> jac(ptr, D, D); ptr += D*D;

          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> invjac(ptr, D, D); ptr += D*D;

          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> xtmp(ptr, D); ptr += D;

          Kokkos::DynRankView<value_type,
                              Kokkos::Impl::ActiveExecutionMemorySpace,
                              Kokkos::MemoryUnmanaged> xold(ptr, D); ptr += D;
   
          // set initial guess
          for (ordinal_type j=0;j<D;++j) xold(j) = 0;
       
          const double tol = tolerence();
          for (ordinal_type iter=0;iter<Parameters::MaxNewton;++iter) {
            // xtmp := F(xold);
            implBasisType::template Serial<OPERATOR_VALUE>::getValues(vals, xold);
            mapToPhysicalFrame(xtmp, vals, nodes);

            // DF^{-1}
            implBasisType::template Serial<OPERATOR_GRAD>::getValues(grads, xold);
            CellTools::Serial::computeJacobian(jac, grads, nodes);
            Kernels::Serial::inverse(invjac, jac);
            
            // Newton
            Kernels::Serial::z_is_axby(xtmp, 1.0, xphy, -1.0, xtmp);  // xtmp := xphy - F(xold);
            Kernels::Serial::gemv_notrans(1.0, invjac, xtmp, 0.0, xref); // xref := DF^{-1}( xphy - F(xold))
            Kernels::Serial::z_is_axby(xref, 1.0, xold,  1.0, xref); // xref += xold
            
            // l2 error
            Kernels::Serial::z_is_axby(xtmp, 1.0, xold, -1.0, xref);

            double err = Kernels::Serial::norm(xtmp, NORM_ONE);

            if (err < tol) 
              break;

            Kernels::Serial::copy(xold, xref);
          }
        }
      };
    };
  }
}

#endif

