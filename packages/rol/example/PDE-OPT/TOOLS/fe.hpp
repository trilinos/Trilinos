// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  fe.hpp
    \brief Given a set of cells with geometric node information, sets up
           data structures used in finite element integration in HGRAD.
*/

#ifndef PDEOPT_FE_HPP
#define PDEOPT_FE_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

template <class Real>
class FE {

private:

  const Teuchos::RCP<Intrepid::FieldContainer<Real> > cellNodes_;                            // coordinates of the cell nodes
  Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basis_;              // Intrepid basis
  const Teuchos::RCP<Intrepid::Cubature<Real, Intrepid::FieldContainer<Real> > > cubature_;  // Intrepid cubature (quadrature, integration) rule
  const int sideId_;                                                                         // local side id for boundary integration

  int c_;    // number of cells in the FE object
  int f_;    // number of basis functions per cell
  int p_;    // number of cubature points per cell
  int d_;    // space dimension of the (parent) cells
  int sd_;   // space dimension of the side cells

  Teuchos::RCP<shards::CellTopology> cellTopo_;   // base (parent) cell topology
  Teuchos::RCP<shards::CellTopology> sideTopo_;   // side (subcell) topology; assumed uniform

  std::vector<std::vector<int> > sideDofs_;       // local dofs on cell sides; 1st index: side number; 2nd index: dof number

  // Containers for local finite element data.
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPoints_;             // points of the cubature rule on the reference cell
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubWeights_;            // weights of the cubature rule on the reference cell
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPointsSubcell_;      // cubature points on the side reference cell
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubWeightsSubcell_;     // cubature weights on the side reference cell
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJac_;               // cell Jacobian matrices
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacInv_;            // inverses of cell Jacobians
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellJacDet_;            // determinants of cell Jacobians
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;   // cell measure (Jacobian determinant) multiplied by the cubature weights
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valReference_;          // value of FE basis in reference space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradReference_;         // gradient of FE basis in reference space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysical_;           // value of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysical_;          // gradient of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalX_;         // x-component of gradient of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalY_;         // y-component of gradient of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalZ_;         // z-component of gradient of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalXWeighted_; // x-component of gradient of FE basis in physical space weighted
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalYWeighted_; // y-component of gradient of FE basis in physical space weighted
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalZWeighted_; // z-component of gradient of FE basis in physical space weighted
  Teuchos::RCP<Intrepid::FieldContainer<Real> > divPhysical_;           // divergence of FE basis in physical space
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;   // value of FE basis in physical space multiplied by weighted cell measure
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradPhysicalWeighted_;  // gradient of FE basis in physical space multiplied by weighted cell measure
  Teuchos::RCP<Intrepid::FieldContainer<Real> > divPhysicalWeighted_;   // divergence of FE basis in physical space multiplied by weighted cell measure
  Teuchos::RCP<Intrepid::FieldContainer<Real> > gradgradMats_;          // cell stiffness matrices
  Teuchos::RCP<Intrepid::FieldContainer<Real> > valvalMats_;            // cell mass matrices
  Teuchos::RCP<Intrepid::FieldContainer<Real> > cubPointsPhysical_;     // cubature points on the physical cells
  Teuchos::RCP<Intrepid::FieldContainer<Real> > dofPoints_;             // degree of freedom points on the reference cell

public:

  FE(const Teuchos::RCP<Intrepid::FieldContainer<Real> >                            & cellNodes,
     const Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > >    & basis,
     const Teuchos::RCP<Intrepid::Cubature<Real, Intrepid::FieldContainer<Real> > > & cubature) :
       cellNodes_(cellNodes), basis_(basis), cubature_(cubature), sideId_(-1) {

    // Get base cell topology from basis.
    cellTopo_ = Teuchos::rcp(new shards::CellTopology(basis_->getBaseCellTopology()));

    // Compute dimensions of multidimensional array members.
    c_  = cellNodes_->dimension(0);
    f_  = basis_->getCardinality();
    p_  = cubature_->getNumPoints();
    d_  = cellTopo_->getDimension();
    sd_ = d_ - 1;

    // Get side subcell topology.
    sideTopo_ = Teuchos::null;

    // Allocate multidimensional arrays.
    cubPoints_               = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_, d_));
    cubWeights_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_));
    cellJac_                 = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_, d_));
    cellJacInv_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_, d_));
    cellJacDet_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_));
    cellWeightedMeasure_     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_));
    valReference_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, p_));
    gradReference_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, p_, d_));
    valPhysical_             = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysical_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_, d_));
    gradPhysicalX_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalY_           = Teuchos::null;
    gradPhysicalZ_           = Teuchos::null;
    if (d_ > 1) {
      gradPhysicalY_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    if (d_ > 2) {
      gradPhysicalZ_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    gradPhysicalXWeighted_   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalYWeighted_   = Teuchos::null;
    gradPhysicalZWeighted_   = Teuchos::null;
    if (d_ > 1) {
      gradPhysicalYWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    if (d_ > 2) {
      gradPhysicalZWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    divPhysical_             = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    valPhysicalWeighted_     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalWeighted_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_, d_));
    divPhysicalWeighted_     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradgradMats_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, f_));
    valvalMats_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, f_));
    cubPointsPhysical_       = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_));
    dofPoints_               = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, d_));

    /*** START: Fill multidimensional arrays. ***/

    // Compute cubature points and weights.
    cubature_->getCubature(*cubPoints_, *cubWeights_);

    // Compute reference basis value and gradient.
    basis_->getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    basis_->getValues(*valReference_, *cubPoints_, Intrepid::OPERATOR_VALUE);       // evaluate value operator at cubature points

    // Compute cell Jacobian matrices, its inverses and determinants.
    Intrepid::CellTools<Real>::setJacobian(*cellJac_, *cubPoints_, *cellNodes_, *cellTopo_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv_, *cellJac_);                       // compute inverses of cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(*cellJacDet_, *cellJac_);                       // compute determinants of cell Jacobians

    // Compute weighted cell measure, i.e., det(J)*(cubature weight).
    Intrepid::FunctionSpaceTools::computeCellMeasure<Real>(*cellWeightedMeasure_,
                                                           *cellJacDet_,
                                                           *cubWeights_);

    // Transform reference values into physical space.
    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(*valPhysical_,
                                                            *valReference_);

    // Multiply with weighted measure to get weighted values in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);

    // Transform reference gradients into physical space.
    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Real>(*gradPhysical_,
                                                           *cellJacInv_,
                                                           *gradReference_);

    // Multiply with weighted measure to get weighted gradients in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*gradPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *gradPhysical_);

    // Extract individual (x, y, z) components of the gradients in physical space.
    for (int c=0; c<c_; ++c) {
      for (int f=0; f<f_; ++f) {
        for (int p=0; p<p_; ++p) {
          (*gradPhysicalX_)(c,f,p) = (*gradPhysical_)(c,f,p,0);
          (*gradPhysicalXWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,0);
          if (d_ > 1) {
          (*gradPhysicalY_)(c,f,p) = (*gradPhysical_)(c,f,p,1);
          (*gradPhysicalYWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,1);
          }
          if (d_ > 2) {
          (*gradPhysicalZ_)(c,f,p) = (*gradPhysical_)(c,f,p,2);
          (*gradPhysicalZWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,2);
          }
        }
      }
    }

    // Build divergence in physical space.
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalX_);
    if (d_ > 1) {
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalY_);
    }
    if (d_ > 2) {
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalZ_);
    }

    // Multiply with weighted measure to get weighted divegence in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*divPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *divPhysical_);

    // Compute stiffness matrices.
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,
                                                  *gradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // Compute mass matrices.
    Intrepid::FunctionSpaceTools::integrate<Real>(*valvalMats_,
                                                  *valPhysical_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // Map reference cubature points to cells in physical space.
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*cubPointsPhysical_,
                                                  *cubPoints_,
                                                  *cellNodes_,
                                                  *cellTopo_);

    // Compute local degrees of freedom on reference cell sides.
    int numSides = cellTopo_->getSideCount();
    if (cellTopo_->getDimension() == 1) {
      numSides = 2;
    }
    if ( numSides ) {
      for (int i=0; i<numSides; ++i) {
        sideDofs_.push_back(computeBoundaryDofs(i));
      }
    }
    else {
      sideDofs_.push_back(computeBoundaryDofs(0));
    }

    // Get coordinates of DOFs in reference cell.
    Teuchos::RCP<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface =
      Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basis_);
    //if (d_ > 1) {
      coord_iface->getDofCoords(*dofPoints_);
    //}

    /*** END: Fill multidimensional arrays. ***/

  }

  FE(const Teuchos::RCP<Intrepid::FieldContainer<Real> >                            & cellNodes,
     const Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > >    & basis,
     const Teuchos::RCP<Intrepid::Cubature<Real, Intrepid::FieldContainer<Real> > > & cubature,
     const int                                                                      & sideId) :
       cellNodes_(cellNodes), basis_(basis), cubature_(cubature), sideId_(sideId) {

    // Get base cell topology from basis.
    cellTopo_ = Teuchos::rcp(new shards::CellTopology(basis_->getBaseCellTopology()));

    // Compute dimensions of multidimensional array members.
    c_  = cellNodes_->dimension(0);
    f_  = basis_->getCardinality();
    p_  = cubature_->getNumPoints();
    d_  = cellTopo_->getDimension();
    sd_ = d_ - 1;
    //std::cout << "FE: c = " << c_ << ", f = " << f_ << ", p = " << p_ << ", d = " << d_ << std::endl;

    // Get side subcell topology.
    sideTopo_ = Teuchos::rcp(new shards::CellTopology(cellTopo_->getCellTopologyData(sd_, sideId_)));

    // Allocate multidimensional arrays.
    cubPoints_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_, d_));
    cubWeights_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_));
    cubPointsSubcell_     = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_, sd_));
    cubWeightsSubcell_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(p_));
    cellJac_              = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_, d_));
    cellJacInv_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_, d_));
    cellJacDet_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_));
    cellWeightedMeasure_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_));
    valReference_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, p_));
    gradReference_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, p_, d_));
    valPhysical_          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysical_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_, d_));
    gradPhysicalX_        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalY_        = Teuchos::null;
    gradPhysicalZ_        = Teuchos::null;
    if (d_ > 1) {
      gradPhysicalY_      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    if (d_ > 2) {
      gradPhysicalZ_      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    divPhysical_          = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    valPhysicalWeighted_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_, d_));
    gradPhysicalXWeighted_   = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradPhysicalYWeighted_   = Teuchos::null;
    gradPhysicalZWeighted_   = Teuchos::null;
    if (d_ > 1) {
      gradPhysicalYWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    if (d_ > 2) {
      gradPhysicalZWeighted_ = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    }
    divPhysicalWeighted_  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, p_));
    gradgradMats_         = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, f_));
    valvalMats_           = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, f_, f_));
    cubPointsPhysical_    = Teuchos::rcp(new Intrepid::FieldContainer<Real>(c_, p_, d_));
    dofPoints_            = Teuchos::rcp(new Intrepid::FieldContainer<Real>(f_, d_));

    /*** START: Fill multidimensional arrays. ***/

    // Compute cubature points and weights.
    cubature_->getCubature(*cubPointsSubcell_, *cubWeights_);

    // Compute reference basis value and gradient.
    Intrepid::CellTools<Real>::mapToReferenceSubcell(*cubPoints_,
                                                     *cubPointsSubcell_,
                                                     sd_,
                                                     sideId_,
                                                     *cellTopo_);
    basis_->getValues(*gradReference_, *cubPoints_, Intrepid::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    basis_->getValues(*valReference_, *cubPoints_, Intrepid::OPERATOR_VALUE);       // evaluate value operator at cubature points

    // Compute cell Jacobian matrices, its inverses and determinants.
    Intrepid::CellTools<Real>::setJacobian(*cellJac_, *cubPoints_, *cellNodes_, *cellTopo_);  // compute cell Jacobians
    Intrepid::CellTools<Real>::setJacobianInv(*cellJacInv_, *cellJac_);                       // compute inverses of cell Jacobians
    Intrepid::CellTools<Real>::setJacobianDet(*cellJacDet_, *cellJac_);                       // compute determinants of cell Jacobians

    // Compute weighted cell measure.
    if (d_ == 2) {
      Intrepid::FunctionSpaceTools::computeEdgeMeasure<Real>(*cellWeightedMeasure_,
                                                             *cellJac_,
                                                             *cubWeights_,
                                                             sideId_,
                                                             *cellTopo_);
    }
    else if (d_ == 3) {
      Intrepid::FunctionSpaceTools::computeFaceMeasure<Real>(*cellWeightedMeasure_,
                                                             *cellJac_,
                                                             *cubWeights_,
                                                             sideId_,
                                                             *cellTopo_);
    }

    // Transform reference values into physical space.
    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(*valPhysical_,
                                                            *valReference_);

    // Multiply with weighted measure to get weighted values in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);

    // Transform reference gradients into physical space.
    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Real>(*gradPhysical_,
                                                           *cellJacInv_,
                                                           *gradReference_);

    // Multiply with weighted measure to get weighted gradients in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*gradPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *gradPhysical_);

    // Extract individual (x, y, z) components of the weighted gradients in physical space.
    for (int c=0; c<c_; ++c) {
      for (int f=0; f<f_; ++f) {
        for (int p=0; p<p_; ++p) {
          (*gradPhysicalX_)(c,f,p) = (*gradPhysical_)(c,f,p,0);
          (*gradPhysicalXWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,0);
          if (d_ > 1) {
          (*gradPhysicalY_)(c,f,p) = (*gradPhysical_)(c,f,p,1);
          (*gradPhysicalYWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,1);
          }
          if (d_ > 2) {
          (*gradPhysicalZ_)(c,f,p) = (*gradPhysical_)(c,f,p,2);
          (*gradPhysicalZWeighted_)(c,f,p) = (*gradPhysicalWeighted_)(c,f,p,2);
          }
        }
      }
    }

    // Build divergence in physical space.
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalX_);
    if (d_ > 1) {
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalY_);
    }
    if (d_ > 2) {
    Intrepid::RealSpaceTools<Real>::add(*divPhysical_, *gradPhysicalZ_);
    }

    // Multiply with weighted measure to get weighted divegence in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*divPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *divPhysical_);

    // Compute stiffness matrices.
    Intrepid::FunctionSpaceTools::integrate<Real>(*gradgradMats_,
                                                  *gradPhysical_,
                                                  *gradPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // Compute mass matrices.
    Intrepid::FunctionSpaceTools::integrate<Real>(*valvalMats_,
                                                  *valPhysical_,
                                                  *valPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // Map reference cubature points to cells in physical space.
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*cubPointsPhysical_,
                                                  *cubPoints_,
                                                  *cellNodes_,
                                                  *cellTopo_);

    // Compute local degrees of freedom on reference cell sides.
    int numSides = cellTopo_->getSideCount();
    for (int i=0; i<numSides; ++i) {
      sideDofs_.push_back(computeBoundaryDofs(i));
    }

    // Get coordinates of DOFs in reference cell.
    Teuchos::RCP<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface =
      Teuchos::rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basis_);
    if (d_ > 1) {
      coord_iface->getDofCoords(*dofPoints_);
    }

    /*** END: Fill multidimensional arrays. ***/

  }

  /** \brief  Returns cell Jacobian matrices at cubature points.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > J() const {
    return cellJac_;
  }

  /** \brief  Returns inverses of cell Jacobians at cubature points.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > invJ() const {
    return cellJacInv_;
  }

  /** \brief  Returns determinants of cell Jacobians at cubature points.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > detJ() const {
    return cellJacDet_;
  }

  /** \brief  Returns values of FE basis at cubature points in reference space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > Nref() const {
    return valReference_;
  }

  /** \brief  Returns gradients of FE basis at cubature points in reference space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > gradNref() const {
    return gradReference_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > N() const {
    return valPhysical_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > NdetJ() const {
    return valPhysicalWeighted_;
  }

  /** \brief  Returns gradient of FE basis at cubature points in physical space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > gradN() const {
    return gradPhysical_;
  }

  /** \brief  Returns gradient of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > gradNdetJ() const {
    return gradPhysicalWeighted_;
  }

  /** \brief  Returns divergence of FE basis at cubature points in physical space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > divN() const {
    return divPhysical_;
  }

  /** \brief  Returns divergence of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > divNdetJ() const {
    return divPhysicalWeighted_;
  }

  /** \brief  Returns x, y or z component of the gradient of FE basis at
              cubature points in physical space.

      \param  coord    [in]   - coordinate index (x=0, y=1, z=2)
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > DND(const int & coord) const {
    if (coord == 0) {
      return gradPhysicalX_;
    }
    else if (coord == 1) {
      return gradPhysicalY_;
    }
    else if (coord == 2) {
      return gradPhysicalZ_;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (PDEOPT::FE::DND): Invalid coordinate argument!");
    }
  }

  /** \brief  Returns x, y or z component of the gradient of FE basis at
              cubature points in physical space, multiplied by weighted
              cell measures.

      \param  coord    [in]   - coordinate index (x=0, y=1, z=2)
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > DNDdetJ(const int & coord) const {
    if (coord == 0) {
      return gradPhysicalXWeighted_;
    }
    else if (coord == 1) {
      return gradPhysicalYWeighted_;
    }
    else if (coord == 2) {
      return gradPhysicalZWeighted_;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (PDEOPT::FE::DNDdetJ): Invalid coordinate argument!");
    }
  }

  /** \brief  Returns stiffness matrices on cells.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > stiffMat() const {
    return gradgradMats_;
  }

  /** \brief  Returns mass matrices on cells.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > massMat() const {
    return valvalMats_;
  }

  /** \brief Returns cubature points on cells in physical space.
  */
  Teuchos::RCP<const Intrepid::FieldContainer<Real> > cubPts() const {
    return cubPointsPhysical_;
  }

  /** \brief Builds FE value interpolant and evaluates it at cubature
             points in physical space.
  */
  void evaluateValue(const Teuchos::RCP<Intrepid::FieldContainer<Real> > & fVals,
                     const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & inCoeffs) const {
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fVals, *inCoeffs, *valPhysical_);
  }

  /** \brief Builds FE gradient interpolant and evaluates it at cubature
             points in physical space.
  */
  void evaluateGradient(const Teuchos::RCP<Intrepid::FieldContainer<Real> > & fGrads,
                        const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & inCoeffs) const {
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fGrads, *inCoeffs, *gradPhysical_);
  }

  /** \brief Computes integral of the product or dot-product of interpolated
             FE fields f1 and f2, indexed by (C,P), for values, or (C,P,D),
             for gradients.
  */
  void computeIntegral(const Teuchos::RCP<Intrepid::FieldContainer<Real> > & integral,
                       const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & f1,
                       const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & f2,
                       const bool sumInto = false) const {
    Intrepid::FieldContainer<Real> f2Weighted(*f2);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(f2Weighted,              // multiply with weighted measure
                                                               *cellWeightedMeasure_,
                                                               *f2);
    Intrepid::FunctionSpaceTools::integrate<Real>(*integral,                            // compute integral of f1*f2
                                                  *f1,
                                                  f2Weighted,
                                                  Intrepid::COMP_CPP,
                                                  sumInto);
  }

  /** \brief Computes the degrees of freedom on a side.
  */
  std::vector<int> computeBoundaryDofs(const int & locSideId) const {

    std::vector<int> bdrydofs;
    std::vector<int> nodeids, edgeids, faceids;

    if (cellTopo_->getDimension() == 1) {
      nodeids.push_back(locSideId);
    }
    else if (cellTopo_->getDimension() == 2) {
      edgeids.push_back(locSideId);
      int numVertices = 2;
      for (int i=0; i<numVertices; ++i) {
        nodeids.push_back(cellTopo_->getNodeMap(1, locSideId, i));
      }
      //for (unsigned i=0; i<nodeids.size(); ++i) {
      //  std::cout << "\nnodeid = " << nodeids[i];
      //}
      //for (unsigned i=0; i<edgeids.size(); ++i) {
      //  std::cout << "\nedgeid = " << edgeids[i];
      //}
    }
    else if (cellTopo_->getDimension() == 3) {
      faceids.push_back(locSideId);
      CellTopologyData_Subcell face = cellTopo_->getCellTopologyData()->subcell[2][locSideId];
      shards::CellTopology face_topo(face.topology);
      int numVertices = face_topo.getVertexCount();
      for (int i=0; i<numVertices; ++i) {
        nodeids.push_back(cellTopo_->getNodeMap(2, locSideId, i));
      }
      int numEdges = face_topo.getEdgeCount();
      for (int i=0; i<numEdges; ++i) {
        edgeids.push_back(mapCellFaceEdge(cellTopo_->getCellTopologyData(), locSideId, i));
      }
      //for (unsigned i=0; i<nodeids.size(); ++i) {
      //  std::cout << "\nnodeid = " << nodeids[i];
      //}
      //for (unsigned i=0; i<edgeids.size(); ++i) {
      //  std::cout << "\nedgeid = " << edgeids[i];
      //}
      //for (unsigned i=0; i<faceids.size(); ++i) {
      //  std::cout << "\nfaceid = " << faceids[i];
      //}
    }

    std::vector<std::vector<std::vector<int> > > tagToId = basis_->getDofOrdinalData();
    std::vector<std::vector<std::vector<int> > > tagToIdCompact;
    int scdim = tagToId.size();
    tagToIdCompact.resize(scdim);
    for (int i=0; i<scdim; ++i) {
      int scid = tagToId[i].size();
      for (int j=0; j<scid; ++j) {
        int doford = tagToId[i][j].size();
        std::vector<int> ids;
        for (int k=0; k<doford; ++k) {
          //std::cout << "\n  i=" << i << "  j=" << j << "  k=" << k << "  id=" << tagToId[i][j][k];
          if (tagToId[i][j][k] != -1) {
            ids.push_back(tagToId[i][j][k]);
          }
        }
        tagToIdCompact[i].push_back(ids);
        int dofordcompact = tagToIdCompact[i][j].size();
        for (int k=0; k<dofordcompact; ++k) {
          //std::cout << "\n  i=" << i << "  j=" << j << "  k=" << k << "  id=" << tagToIdCompact[i][j][k];
        }
      }
    }

    int numNodeIds = nodeids.size();
    if (tagToIdCompact.size() > 0) {
      for (int i=0; i<numNodeIds; ++i) {
        int numdofs = tagToIdCompact[0][nodeids[i]].size();
        for (int j=0; j<numdofs; ++j) {
          bdrydofs.push_back(tagToIdCompact[0][nodeids[i]][j]);
        }
      }
    }

    int numEdgeIds = edgeids.size();
    if (tagToIdCompact.size() > 1) {
      for (int i=0; i<numEdgeIds; ++i) {
        int numdofs = tagToIdCompact[1][edgeids[i]].size();
        for (int j=0; j<numdofs; ++j) {
          bdrydofs.push_back(tagToIdCompact[1][edgeids[i]][j]);
        }
      }
    }

    int numFaceIds = faceids.size();
    if (tagToIdCompact.size() > 2) {
      for (int i=0; i<numFaceIds; ++i) {
        int numdofs = tagToIdCompact[2][faceids[i]].size();
        for (int j=0; j<numdofs; ++j) {
          bdrydofs.push_back(tagToIdCompact[2][faceids[i]][j]);
        }
      }
    }

    //for (unsigned i=0; i<bdrydofs.size(); ++i) {
    //  std::cout << "\ndofid = " << bdrydofs[i];
    //}

    return bdrydofs;

  }

  /** \brief Returns the degrees of freedom on reference cell sides.
  */
  const std::vector<std::vector<int> > & getBoundaryDofs() const {
    return sideDofs_;
  }

  /** \brief Computes coordinates of degrees of freedom on cells.
  */
  void computeDofCoords(const Teuchos::RCP<Intrepid::FieldContainer<Real> > & dofCoords,
                        const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & cellNodes) const {
    // Map reference DOF locations to physical space.
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*dofCoords,
                                                  *dofPoints_,
                                                  *cellNodes,
                                                  *cellTopo_);
  }

}; // FE

#endif
