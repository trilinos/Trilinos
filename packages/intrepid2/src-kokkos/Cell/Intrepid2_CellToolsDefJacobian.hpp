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


/** \file   Intrepid_CellToolsDef.hpp
    \brief  Definition file for the Intrepid2::CellTools class.
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

  template<typename SpT>
  template<typename jacobianValueType, class ...jacobianProperties,
           typename pointValueType,    class ...pointProperties,
           typename worksetValueType,  class ...worksetProperties>
  void
  CellTools<SpT>::
  setJacobian( /**/  Kokkos::DynRankView<jacobianValueType,jacobianProperties...> jacobian,
               const Kokkos::DynRankView<pointValueType,pointProperties...>       points,
               const Kokkos::DynRankView<worksetValueType,worksetProperties...>   workset,
               const shards::CellTopology cellTopo ) {
#ifdef HAVE_INTREPID2_DEBUG    
    validateArguments_setJacobian(jacobian, points, cellWorkset, whichCell,  cellTopo);
#endif
    
    const auto spaceDim = cellTopo.getDimension();
    const auto numCells = workset.dimension(0);
    
    //points can be rank-2 (P,D), or rank-3 (C,P,D)
    const auto pointRank = points.rank();
    const auto numPoints = (pointRank == 2 ? points.dimension(0) : points.dimension(1));
    
    // Jacobian is computed using gradients of an appropriate H(grad) basis function
    const auto basis = getHgradBasis(cellTopo);
    const auto basisCardinality = basis->getCardinality();    

    typedef typename ExecSpace<typename whatever::execution_space,SpT>::ExecSpaceType ExecSpaceType;
    
    // Temp (F,P,D) array for the values of basis functions gradients at the reference points
    // this needs to be thread private; range policy does not have this capability yet...
    // use manual ... for now.... we do not have team interface...
    // let's use team interface.
    //Kokkos::DynRankView<scratchValueType,scratchProperties...> gradients(scratch.data(), 
    //                                                                     basisCardinality,
    //                                                                     numPoints,
    //                                                                     spaceDim);

    //(*basis).Serial<OPERATOR_GRAD>::getValues(gradients, );

    const auto loopSize = numCells;
    Kokkos::TeamPolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(jacobian, ) );
    
    
    //     // refPoints is (P,D): a single or multiple cell jacobians computed for a single set of ref. points
    //   (*basis).getValues(gradients, points, OPERATOR_GRAD);
    //   break;
    // }
    //   case 3:
    //     // points is (C,P,D): multiple jacobians computed at multiple point sets, one jacobian per cell  

    //     break;
    //   }  
    // default:
    //       INTREPID2_TEST_FOR_EXCEPTION( !( (getrank(points) == 2) && (getrank(points) == 3) ), std::invalid_argument,
    //                           ">>> ERROR (Intrepid2::CellTools::setJacobian): rank 2 or 3 required for points array. ");        
    //   }//switch
  }

  template<typename SpT>
  template<typename jacobianInvValueType, class ...jacobianInvProperties,                                   
           typename jacobianValueType,    class ...jacobianProperties>                                      
  void                                                                                               
  CellTools<SpT>::
  setJacobianInv( /**/  Kokkos::DynRankView<jacobianInvValueType,jacobianInvProperties...> jacobianInv,     
                  const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>       jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    validateArguments_setJacobianInv(jacobianInv, jacobian);
#endif
    RealSpaceTools<SpT>::inverse(jacobianInv, jacobian);
  }
  
  template<typename SpT>
  template<typename jacobianDetValueType, class ...jacobianDetProperties,                                   
           typename jacobianValueType,    class ...jacobianProperties>                                      
  void                                                                                               
  CellTools<SpT>::
  setJacobianDet( /**/  Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...>  jacobianDet,    
                  const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>        jacobian ) {
#ifdef HAVE_INTREPID2_DEBUG
    validateArguments_setJacobianDetArgs(jacobianDet, jacobian);
#endif
    RealSpaceTools<SpT>::det(jacobianDet, jacobian);
  }

}

#endif
