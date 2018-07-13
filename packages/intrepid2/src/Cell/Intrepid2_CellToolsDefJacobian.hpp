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


/** \file   Intrepid2_CellToolsDefJacobian.hpp
    \brief  Definition file for the Jacobian functions in the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/
#ifndef __INTREPID2_CELLTOOLS_DEF_HPP__
#define __INTREPID2_CELLTOOLS_DEF_HPP__

// disable clang warnings
#if defined (__clang__) && !defined (__INTEL_COMPILER)
#pragma clang system_header
#endif

namespace Intrepid2 {


  //============================================================================================//
  //                                                                                            //
  //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
  //                                                                                            //
  //============================================================================================//

  namespace FunctorCellTools {
    /**
     \brief Functor for calculation of Jacobian on cell workset see Intrepid2::CellTools for more
    */
    template<typename jacobianViewType, 
             typename worksetCellType, 
             typename basisGradType>
    struct F_setJacobian {
            jacobianViewType _jacobian;
      const worksetCellType  _worksetCells;
      const basisGradType    _basisGrads;

      KOKKOS_INLINE_FUNCTION
      F_setJacobian( jacobianViewType jacobian_,
                     worksetCellType  worksetCells_,
                     basisGradType    basisGrads_ )
        : _jacobian(jacobian_), _worksetCells(worksetCells_), _basisGrads(basisGrads_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl, 
                      const ordinal_type pt) const {
              auto jac  = Kokkos::subview( _jacobian,     cl, pt, Kokkos::ALL(), Kokkos::ALL());
        const auto dofs = Kokkos::subview( _worksetCells, cl,     Kokkos::ALL(), Kokkos::ALL());
        
        const ordinal_type gradRank = _basisGrads.rank(); 
        const auto grad = ( gradRank == 3 ? Kokkos::subdynrankview( _basisGrads,     Kokkos::ALL(), pt, Kokkos::ALL()) :
                                            Kokkos::subdynrankview( _basisGrads, cl, Kokkos::ALL(), pt, Kokkos::ALL()));
        
        const ordinal_type dim = jac.extent(0); // dim0 and dim1 should match
        const ordinal_type cardinality = grad.extent(0);

        for (ordinal_type i=0;i<dim;++i)
          for (ordinal_type j=0;j<dim;++j) {
            jac(i, j) = 0;
            for (ordinal_type bf=0;bf<cardinality;++bf) 
              jac(i, j) += dofs(bf, i)*grad(bf, j);
          } 
      }
    };
  }

  template<typename SpT>
  template<typename jacobianValueType,    class ...jacobianProperties,
           typename pointValueType,       class ...pointProperties,
           typename worksetCellValueType, class ...worksetCellProperties,
           typename HGradBasisPtrType>
  void
  CellTools<SpT>::
  setJacobian(       Kokkos::DynRankView<jacobianValueType,jacobianProperties...>       jacobian,
               const Kokkos::DynRankView<pointValueType,pointProperties...>             points,
               const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
               const HGradBasisPtrType basis ) {
#ifdef HAVE_INTREPID2_DEBUG    
    CellTools_setJacobianArgs(jacobian, points, worksetCell, basis->getBaseCellTopology());
    //static_assert(std::is_same( pointValueType, decltype(basis->getDummyOutputValue()) ));
#endif
    const auto cellTopo = basis->getBaseCellTopology();
    const ordinal_type spaceDim = cellTopo.getDimension();
    const ordinal_type numCells = worksetCell.extent(0);
    
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    const ordinal_type pointRank = points.rank();
    const ordinal_type numPoints = (pointRank == 2 ? points.extent(0) : points.extent(1));
    const ordinal_type basisCardinality = basis->getCardinality();    

    typedef Kokkos::DynRankView<jacobianValueType,jacobianProperties...> jacobianViewType;
    
    // the following does not work for RCP; its * operator returns reference not the object
    //typedef typename decltype(*basis)::output_value_type gradValueType;
    //typedef Kokkos::DynRankView<decltype(basis->getDummyOutputValue()),SpT> gradViewType;

    auto vcprop = Kokkos::common_view_alloc_prop(points);
    typedef Kokkos::DynRankView<typename decltype(vcprop)::value_type,SpT> gradViewType;

    gradViewType grads;

    typedef Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCellViewType;
    typedef typename ExecSpace<typename worksetCellViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    switch (pointRank) {
    case 2: {
      // For most FEMs
      grads = gradViewType(Kokkos::view_alloc("CellTools::setJacobian::grads", vcprop),basisCardinality, numPoints, spaceDim);
      basis->getValues(grads, 
                       points, 
                       OPERATOR_GRAD);
      break;
    }
    case 3: { 
      // For CVFEM
      grads = gradViewType(Kokkos::view_alloc("CellTools::setJacobian::grads", vcprop), numCells, basisCardinality, numPoints, spaceDim);
      for (ordinal_type cell=0;cell<numCells;++cell) 
        basis->getValues(Kokkos::subview( grads,  cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL() ),  
                         Kokkos::subview( points, cell, Kokkos::ALL(), Kokkos::ALL() ),  
                         OPERATOR_GRAD);
      break;
    }
    }

    typedef Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCellViewType;
    typedef FunctorCellTools::F_setJacobian<jacobianViewType,worksetCellViewType,gradViewType> FunctorType;
    typedef typename ExecSpace<typename worksetCellViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;
    
    using range_policy_type = Kokkos::Experimental::MDRangePolicy
      < ExecSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    range_policy_type policy( { 0, 0 },
                              { jacobian.extent(0), jacobian.extent(1) } );
    Kokkos::parallel_for( policy, FunctorType(jacobian, worksetCell, grads) );
  }

  template<typename SpT>
  template<typename jacobianInvValueType, class ...jacobianInvProperties,                                   
           typename jacobianValueType,    class ...jacobianProperties>                                      
  void                                                                                               
  CellTools<SpT>::
  setJacobianInv(       Kokkos::DynRankView<jacobianInvValueType,jacobianInvProperties...> jacobianInv,     
                  const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>       jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_setJacobianInvArgs(jacobianInv, jacobian);
#endif
    RealSpaceTools<SpT>::inverse(jacobianInv, jacobian);
  }
  
  template<typename SpT>
  template<typename jacobianDetValueType, class ...jacobianDetProperties,                                   
           typename jacobianValueType,    class ...jacobianProperties>                                      
  void                                                                                               
  CellTools<SpT>::
  setJacobianDet(       Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...>  jacobianDet,    
                  const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>        jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    CellTools_setJacobianDetArgs(jacobianDet, jacobian);
#endif
    RealSpaceTools<SpT>::det(jacobianDet, jacobian);
  }

}

#endif
