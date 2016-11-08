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

/** \file test_01.cpp
\brief  Performance test comparing dynrankview overhead
\author Created by Kyungjoo Kim.
*/

#include "Intrepid2_config.h"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"


#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_DefaultCubatureFactory.hpp"

#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"


#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Intrepid2 {
  
  namespace Test {

    template<typename ValueType,
             typename DeviceSpaceType>
    struct F_hgrad_eval {
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> ViewType;
      typedef typename ViewType::const_type constViewType;
      typedef ValueType value_type;

      /**/  ViewType _weighted_basis_values;
      /**/  ViewType _weighted_basis_grads;
      const ViewType _grads;
      const ViewType _workset;
      const ViewType _weights;
      const ViewType _basis_values;
      const ViewType _basis_grads;

      KOKKOS_INLINE_FUNCTION
      F_hgrad_eval(ViewType weighted_basis_values_,
                   ViewType weighted_basis_grads_,
                   const ViewType grads_, // gradient for jacobian
                   const ViewType workset_, // workset input
                   const ViewType weights_, // weights
                   const ViewType basis_values_, // reference values
                   const ViewType basis_grads_) // reference grads
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
        value_type buf[9 + 9];

        const auto grad = Kokkos::subview(_grads,       Kokkos::ALL(), pt, Kokkos::ALL());
        const auto dofs = Kokkos::subview(_workset, cl, Kokkos::ALL(),     Kokkos::ALL());

        const ordinal_type card = dofs.dimension(0);
        const ordinal_type dim = dofs.dimension(1);
        
        // temporary values
        Kokkos::View<value_type**, Kokkos::Impl::ActiveExecutionMemorySpace> 
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
        const value_type det = ( jac(0,0) * jac(1,1) -
                                 jac(0,1) * jac(1,0) );
        
        {
          const value_type val = det;
          jac_inv(0,0) =   jac(1,1)/val;
          jac_inv(1,1) =   jac(0,0)/val;
          
          jac_inv(1,0) = - jac(1,0)/val;
          jac_inv(0,1) = - jac(0,1)/val;
        }
        
        // computeCellMeasure
        const value_type cell_measure = (det > 0 ? det : -det)*_weights(pt);
        
        // multiplyMeasure
        for (ordinal_type bf=0;bf<card;++bf)        
          _weighted_basis_values(cl, bf, pt) = cell_measure*_basis_values(bf, pt);
        
        // // HGRADtransformGRAD
        // auto weighted_grad = Kokkos::subview(_weighted_basis_grads, cl, Kokkos::ALL(), pt, Kokkos::ALL());
        // const auto basis_grad = Kokkos::subview(_basis_grads, Kokkos::ALL(), pt, Kokkos::ALL());
        
        // {
        //   const ordinal_type card = basis_grad.dimension(0);
        //   const ordinal_type dim = basis_grad.dimension(1);
        //   for (ordinal_type bf=0;bf<card;++bf)
        //     for (ordinal_type i=0;i<dim;++i) {
        //       weighted_grad(bf, i) = 0;
        //       for (ordinal_type j=0;j<dim;++j) 
        //         weighted_grad(bf, i) += jac_inv(i,j)*basis_grad(bf, j);
        //       weighted_grad(bf, i) *= cell_measure;
        //     }
        // }
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl, 
                      const ordinal_type pt) const {
        apply(cl, pt);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl) const {
        const ordinal_type npts = _basis_values.dimension(1);
        for (ordinal_type pt=0;pt<npts;++pt) 
          apply(cl, pt);
      }
    };
    
    template<typename ValueType, typename DeviceSpaceType>
    int ComputeBasis_HGRAD(const ordinal_type nworkset,
                           const ordinal_type C,
                           const ordinal_type order,
                           const bool verbose) {

      Teuchos::RCP<std::ostream> verboseStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        verboseStream = Teuchos::rcp(&std::cout, false);
      else
        verboseStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *verboseStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*verboseStream, false);
      *verboseStream << "HostSpace::    ";   HostSpaceType::print_configuration(*verboseStream, false);

      using BasisTypeHost = Basis_HGRAD_QUAD_C1_FEM<HostSpaceType,ValueType,ValueType>;
      using ImplBasisType = Impl::Basis_HGRAD_QUAD_C1_FEM;
      using range_type = Kokkos::pair<ordinal_type,ordinal_type>;

      
      Kokkos::Impl::Timer timer;
      double t_horizontal = 0, t_vertical = 0, t_vectorize = 0;
      int errorFlag = 0;

      BasisTypeHost hostBasis;
      const auto cellTopo = hostBasis.getBaseCellTopology();

      auto cubature = DefaultCubatureFactory::create<DeviceSpaceType,ValueType,ValueType>(cellTopo, order);
      const ordinal_type 
        numCells = C, 
        numVerts = cellTopo.getVertexCount(),
        numDofs = hostBasis.getCardinality(),
        numPoints = cubature->getNumPoints(), 
        spaceDim = cubature->getDimension();

      Kokkos::DynRankView<ValueType,HostSpaceType> dofCoordsHost("dofCoordsHost", numDofs, spaceDim);
      hostBasis.getDofCoords(dofCoordsHost);
      const auto refNodesHost = Kokkos::subview(dofCoordsHost, range_type(0, numVerts), Kokkos::ALL());
      
      // pertub nodes
      Kokkos::DynRankView<ValueType,HostSpaceType> worksetCellsHost("worksetCellsHost", numCells, numVerts, spaceDim);
      for (ordinal_type cell=0;cell<numCells;++cell) {
        for (ordinal_type i=0;i<numVerts;++i)
          for (ordinal_type j=0;j<spaceDim;++j) {
            ValueType val = (rand()/(RAND_MAX + 1.0))*0.2 -0.1;
            worksetCellsHost(cell, i, j) = refNodesHost(i, j) + val; 
          }
      }

      auto worksetCells = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), worksetCellsHost);
      Kokkos::deep_copy(worksetCells, worksetCellsHost);

      Kokkos::DynRankView<ValueType,DeviceSpaceType> refPoints("refPoints", numPoints, spaceDim), refWeights("refWeights", numPoints);
      cubature->getCubature(refPoints, refWeights);

      std::cout
        << "===============================================================================\n" 
        << " Performance Test evaluating ComputeBasis \n"
        << " # of workset = " << nworkset << "\n" 
        << " Test Array Structure (C,F,P,D) = " << numCells << ", " << numDofs << ", " << numPoints << ", " << spaceDim << "\n"
        << "===============================================================================\n";

      *verboseStream
        << "\n"
        << "===============================================================================\n"
        << "TEST 1: evaluateFields \n"
        << "===============================================================================\n";
      
      try {

        
        Kokkos::DynRankView<ValueType,DeviceSpaceType> 
          refBasisValues("refBasisValues", numDofs, numPoints),
          refBasisGrads ("refBasisGrads",  numDofs, numPoints, spaceDim);
        
        ImplBasisType::getValues<DeviceSpaceType>(refBasisValues, refPoints, OPERATOR_VALUE);
        ImplBasisType::getValues<DeviceSpaceType>(refBasisGrads,  refPoints, OPERATOR_GRAD);
        
        const ordinal_type ibegin = -3;
        // testing sequential appraoch
        {          
          Kokkos::DynRankView<ValueType,DeviceSpaceType> 
            jacobian           ("jacobian",            numCells,          numPoints, spaceDim, spaceDim),
            jacobianInv        ("jacobianInv",         numCells,          numPoints, spaceDim, spaceDim),
            jacobianDet        ("jacobianDet",         numCells,          numPoints),
            cellMeasure        ("cellMeasure",         numCells,          numPoints),
            phyBasisValues     ("phyBasisValues",      numCells, numDofs, numPoints),
            phyBasisGrads      ("phyBasisGrads",       numCells, numDofs, numPoints, spaceDim),
            weightedBasisValues("weightedBasisValues", numCells, numDofs, numPoints),
            weightedBasisGrads ("weightedBasisGrads",  numCells, numDofs, numPoints, spaceDim);
          
          typedef CellTools<DeviceSpaceType> cts;
          typedef FunctionSpaceTools<DeviceSpaceType> fts;


          for (ordinal_type iwork=ibegin;iwork<nworkset;++iwork) {
            DeviceSpaceType::fence();
            timer.reset();

            cts::setJacobian(jacobian, refPoints, worksetCells, cellTopo);
            cts::setJacobianInv(jacobianInv, jacobian);
            cts::setJacobianDet(jacobianDet, jacobian);
            fts::computeCellMeasure(cellMeasure, jacobianDet, refWeights);
            
            fts::HGRADtransformVALUE(phyBasisValues, refBasisValues);
            fts::multiplyMeasure(weightedBasisValues, cellMeasure, phyBasisValues);
            
            fts::HGRADtransformGRAD(phyBasisGrads, jacobianInv, refBasisGrads);
            fts::multiplyMeasure(weightedBasisGrads, cellMeasure, phyBasisGrads);

            DeviceSpaceType::fence();
            t_horizontal += (iwork >= 0)*timer.seconds();
          }

        }

        // testing vertical approach
        {          
          Kokkos::DynRankView<ValueType,DeviceSpaceType> 
            weightedBasisValues("weightedBasisValues", numCells, numDofs, numPoints),
            weightedBasisGrads ("weightedBasisGrads",  numCells, numDofs, numPoints, spaceDim);

          typedef CellTools<DeviceSpaceType> cts;
          typedef FunctionSpaceTools<DeviceSpaceType> fts;

          typedef F_hgrad_eval<ValueType,DeviceSpaceType> FunctorType;

          //Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, numCells);
          using range_policy_type = Kokkos::Experimental::MDRangePolicy
            < DeviceSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
          range_policy_type policy( {        0,         0 },
                                    { numCells, numPoints } );

          FunctorType functor(weightedBasisValues,
                              weightedBasisGrads,
                              refBasisGrads,
                              worksetCells,
                              refWeights,
                              refBasisValues,
                              refBasisGrads);
          
          for (ordinal_type iwork=ibegin;iwork<nworkset;++iwork) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::Experimental::md_parallel_for(policy, functor);

            DeviceSpaceType::fence();
            t_vertical += (iwork >= 0)*timer.seconds();
          }

        }

      } catch (std::exception err) {
        *verboseStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *verboseStream << err.what() << '\n';
        *verboseStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      std::cout 
        << "TEST HGRAD " 
        << ": t_horizontal = " << (t_horizontal/nworkset) 
        << ", t_vertical = " << (t_vertical/nworkset)
        << ", t_vectorize = " << (t_vectorize/nworkset)
        << std::endl;
      
      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  } // end of namespace TEST
} // end of namespace Intrepid2
