// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  feK.hpp
    \brief Given a set of cells with geometric node information, sets up
           data structures used in finite element integration in HGRAD.
*/

#ifndef PDEOPT_FE_HPP
#define PDEOPT_FE_HPP

#include "ROL_Ptr.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Intrepid2_Basis.hpp"
#include "Intrepid2_Cubature.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

template <class Real, class DeviceType>
class FE {

public:
  using scalar_view = Kokkos::DynRankView<Real,DeviceType>;
  using fst = Intrepid2::FunctionSpaceTools<DeviceType>;
  using ct = Intrepid2::CellTools<DeviceType>;
  using rst = Intrepid2::RealSpaceTools<DeviceType>;

private:

  const scalar_view cellNodes_;                                 // coordinates of the cell nodes
  Intrepid2::BasisPtr<DeviceType, Real, Real> basis_;           // Intrepid basis
  const ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> cubature_;  // Intrepid cubature (quadrature, integration) rule
  const int sideId_;                                                                         // local side id for boundary integration

  int c_;    // number of cells in the FE object
  int f_;    // number of basis functions per cell
  int p_;    // number of cubature points per cell
  int d_;    // space dimension of the (parent) cells
  int sd_;   // space dimension of the side cells

  ROL::Ptr<shards::CellTopology> cellTopo_;   // base (parent) cell topology
  ROL::Ptr<shards::CellTopology> sideTopo_;   // side (subcell) topology; assumed uniform

  std::vector<std::vector<int> > sideDofs_;       // local dofs on cell sides; 1st index: side number; 2nd index: dof number

  // Containers for local finite element data.
  scalar_view cubPoints_;             // points of the cubature rule on the reference cell
  scalar_view cubWeights_;            // weights of the cubature rule on the reference cell
  scalar_view cubPointsSubcell_;      // cubature points on the side reference cell
  scalar_view cubWeightsSubcell_;     // cubature weights on the side reference cell
  scalar_view cellJac_;               // cell Jacobian matrices
  scalar_view cellJacInv_;            // inverses of cell Jacobians
  scalar_view cellJacDet_;            // determinants of cell Jacobians
  scalar_view cellWeightedMeasure_;   // cell measure (Jacobian determinant) multiplied by the cubature weights
  scalar_view valReference_;          // value of FE basis in reference space
  scalar_view gradReference_;         // gradient of FE basis in reference space
  scalar_view valPhysical_;           // value of FE basis in physical space
  scalar_view gradPhysical_;          // gradient of FE basis in physical space
  scalar_view gradPhysicalX_;         // x-component of gradient of FE basis in physical space
  scalar_view gradPhysicalY_;         // y-component of gradient of FE basis in physical space
  scalar_view gradPhysicalZ_;         // z-component of gradient of FE basis in physical space
  scalar_view gradPhysicalXWeighted_; // x-component of gradient of FE basis in physical space weighted
  scalar_view gradPhysicalYWeighted_; // y-component of gradient of FE basis in physical space weighted
  scalar_view gradPhysicalZWeighted_; // z-component of gradient of FE basis in physical space weighted
  scalar_view divPhysical_;           // divergence of FE basis in physical space
  scalar_view valPhysicalWeighted_;   // value of FE basis in physical space multiplied by weighted cell measure
  scalar_view gradPhysicalWeighted_;  // gradient of FE basis in physical space multiplied by weighted cell measure
  scalar_view divPhysicalWeighted_;   // divergence of FE basis in physical space multiplied by weighted cell measure
  scalar_view gradgradMats_;          // cell stiffness matrices
  scalar_view valvalMats_;            // cell mass matrices
  scalar_view cubPointsPhysical_;     // cubature points on the physical cells
  scalar_view dofPoints_;             // degree of freedom points on the reference cell

public:

  FE(const scalar_view                                           cellNodes,
     const Intrepid2::BasisPtr<DeviceType, Real, Real>           basis,
     const ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> cubature,
     bool computeBdryDofs = true) :
       cellNodes_(cellNodes), basis_(basis), cubature_(cubature), sideId_(-1) {

    // Get base cell topology from basis.
    cellTopo_ = ROL::makePtr<shards::CellTopology>(basis_->getBaseCellTopology());

    // Compute dimensions of multidimensional array members.
    c_  = cellNodes_.extent_int(0);
    f_  = basis_->getCardinality();
    p_  = cubature_->getNumPoints();
    d_  = cellTopo_->getDimension();
    sd_ = d_ - 1;

    // Get side subcell topology.
    sideTopo_ = ROL::nullPtr;

    // Allocate multidimensional arrays.
    cubPoints_               = scalar_view("cubPoints_", p_, d_);
    cubWeights_              = scalar_view("cubWeights_", p_);
    cellJac_                 = scalar_view("cellJac_", c_, p_, d_, d_);
    cellJacInv_              = scalar_view("cellJacInv_", c_, p_, d_, d_);
    cellJacDet_              = scalar_view("cellJacDet_", c_, p_);
    cellWeightedMeasure_     = scalar_view("cellWeightedMeasure_", c_, p_);
    valReference_            = scalar_view("valReference_", f_, p_);
    gradReference_           = scalar_view("gradReference_", f_, p_, d_);
    valPhysical_             = scalar_view("valPhysical_", c_, f_, p_);
    gradPhysical_            = scalar_view("gradPhysical_", c_, f_, p_, d_);
    gradPhysicalX_           = scalar_view("gradPhysicalX_", c_, f_, p_);
    if (d_ > 1) {
      gradPhysicalY_         = scalar_view("gradPhysicalY_", c_, f_, p_);
    }
    if (d_ > 2) {
      gradPhysicalZ_         = scalar_view("gradPhysicalZ_", c_, f_, p_);
    }
    gradPhysicalXWeighted_   = scalar_view("gradPhysicalXWeighted_", c_, f_, p_);
    if (d_ > 1) {
      gradPhysicalYWeighted_ = scalar_view("gradPhysicalYWeighted_", c_, f_, p_);
    }
    if (d_ > 2) {
      gradPhysicalZWeighted_ = scalar_view("gradPhysicalZWeighted_", c_, f_, p_);
    }
    divPhysical_             = scalar_view("divPhysical_", c_, f_, p_);
    valPhysicalWeighted_     = scalar_view("valPhysicalWeighted_", c_, f_, p_);
    gradPhysicalWeighted_    = scalar_view("gradPhysicalWeighted_", c_, f_, p_, d_);
    divPhysicalWeighted_     = scalar_view("divPhysicalWeighted_", c_, f_, p_);
    gradgradMats_            = scalar_view("gradgradMats_", c_, f_, f_);
    valvalMats_              = scalar_view("valvalMats_", c_, f_, f_);
    cubPointsPhysical_       = scalar_view("cubPointsPhysical_", c_, p_, d_);
    dofPoints_               = scalar_view("dofPoints_", f_, d_);

    /*** START: Fill multidimensional arrays. ***/

    // Compute cubature points and weights.
    cubature_->getCubature(cubPoints_, cubWeights_);

    // Compute reference basis value and gradient.
    basis_->getValues(gradReference_, cubPoints_, Intrepid2::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    basis_->getValues(valReference_, cubPoints_, Intrepid2::OPERATOR_VALUE);       // evaluate value operator at cubature points

    // Compute cell Jacobian matrices, its inverses and determinants.
    ct::setJacobian(cellJac_, cubPoints_, cellNodes_, *cellTopo_);  // compute cell Jacobians
    ct::setJacobianInv(cellJacInv_, cellJac_);                       // compute inverses of cell Jacobians
    ct::setJacobianDet(cellJacDet_, cellJac_);                       // compute determinants of cell Jacobians

    // Compute weighted cell measure, i.e., det(J)*(cubature weight).
    fst::computeCellMeasure(cellWeightedMeasure_, cellJacDet_, cubWeights_);

    // Transform reference values into physical space.
    fst::HGRADtransformVALUE(valPhysical_, valReference_);

    // Multiply with weighted measure to get weighted values in physical space.
    fst::multiplyMeasure(valPhysicalWeighted_, cellWeightedMeasure_, valPhysical_);

    // Transform reference gradients into physical space.
    fst::HGRADtransformGRAD(gradPhysical_, cellJacInv_, gradReference_);

    // Multiply with weighted measure to get weighted gradients in physical space.
    fst::multiplyMeasure(gradPhysicalWeighted_, cellWeightedMeasure_, gradPhysical_);

    // Extract individual (x, y, z) components of the gradients in physical space.
    for (int c=0; c<c_; ++c) {
      for (int f=0; f<f_; ++f) {
        for (int p=0; p<p_; ++p) {
          gradPhysicalX_(c,f,p) = gradPhysical_(c,f,p,0);
            gradPhysicalXWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,0);
          if (d_ > 1) {
            gradPhysicalY_(c,f,p) = gradPhysical_(c,f,p,1);
            gradPhysicalYWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,1);
          }
          if (d_ > 2) {
            gradPhysicalZ_(c,f,p) = gradPhysical_(c,f,p,2);
            gradPhysicalZWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,2);
          }
        }
      }
    }

    // Build divergence in physical space.
    rst::add(divPhysical_, gradPhysicalX_);
    if (d_ > 1) {
      rst::add(divPhysical_, gradPhysicalY_);
    }
    if (d_ > 2) {
      rst::add(divPhysical_, gradPhysicalZ_);
    }

    // Multiply with weighted measure to get weighted divegence in physical space.
    fst::multiplyMeasure(divPhysicalWeighted_, cellWeightedMeasure_, divPhysical_);

    // Compute stiffness matrices.
    fst::integrate(gradgradMats_, gradPhysical_, gradPhysicalWeighted_);

    // Compute mass matrices.
    fst::integrate(valvalMats_, valPhysical_, valPhysicalWeighted_);

    // Map reference cubature points to cells in physical space.
    ct::mapToPhysicalFrame(cubPointsPhysical_, cubPoints_, cellNodes_, *cellTopo_);

    // Compute local degrees of freedom on reference cell sides.
    if (computeBdryDofs) {
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
    }

    // Get coordinates of DOFs in reference cell.
    basis->getDofCoords(dofPoints_);

    /*** END: Fill multidimensional arrays. ***/

  }

  FE(const scalar_view cellNodes,
     const Intrepid2::BasisPtr<DeviceType, Real, Real>           basis,
     const ROL::Ptr<Intrepid2::Cubature<DeviceType, Real, Real>> cubature,
     const int&                                                  sideId) :
       cellNodes_(cellNodes), basis_(basis), cubature_(cubature), sideId_(sideId) {

    // Get base cell topology from basis.
    cellTopo_ = ROL::makePtr<shards::CellTopology>(basis_->getBaseCellTopology());

    // Compute dimensions of multidimensional array members.
    c_  = cellNodes_.extent_int(0);
    f_  = basis_->getCardinality();
    p_  = cubature_->getNumPoints();
    d_  = cellTopo_->getDimension();
    sd_ = d_ - 1;
    //std::cout << "FE: c = " << c_ << ", f = " << f_ << ", p = " << p_ << ", d = " << d_ << std::endl;

    // Get side subcell topology.
    sideTopo_ = ROL::makePtr<shards::CellTopology>(cellTopo_->getCellTopologyData(sd_, sideId_));

    // Allocate multidimensional arrays.
    cubPoints_            = scalar_view("cubPoints_", p_, d_);
    cubWeights_           = scalar_view("cubWeights_", p_);
    cubPointsSubcell_     = scalar_view("cubPointsSubcell_", p_, sd_);
    cubWeightsSubcell_    = scalar_view("cubWeightsSubcell_", p_);
    cellJac_              = scalar_view("cellJac_", c_, p_, d_, d_);
    cellJacInv_           = scalar_view("cellJacInv_", c_, p_, d_, d_);
    cellJacDet_           = scalar_view("cellJacDet_", c_, p_);
    cellWeightedMeasure_  = scalar_view("cellWeightedMeasure_", c_, p_);
    valReference_         = scalar_view("valReference_", f_, p_);
    gradReference_        = scalar_view("gradReference_", f_, p_, d_);
    valPhysical_          = scalar_view("valPhysical_", c_, f_, p_);
    gradPhysical_         = scalar_view("gradPhysical_", c_, f_, p_, d_);
    gradPhysicalX_        = scalar_view("gradPhysicalX_", c_, f_, p_);
    if (d_ > 1) {
      gradPhysicalY_      = scalar_view("gradPhysicalY_", c_, f_, p_);
    }
    if (d_ > 2) {
      gradPhysicalZ_      = scalar_view("gradPhysicalZ_", c_, f_, p_);
    }
    divPhysical_          = scalar_view("divPhysical_", c_, f_, p_);
    valPhysicalWeighted_  = scalar_view("valPhysicalWeighted_", c_, f_, p_);
    gradPhysicalWeighted_ = scalar_view("gradPhysicalWeighted_", c_, f_, p_, d_);
    gradPhysicalXWeighted_  = scalar_view("gradPhysicalXWeighted_", c_, f_, p_);
    if (d_ > 1) {
      gradPhysicalYWeighted_ = scalar_view("gradPhysicalYWeighted_", c_, f_, p_);
    }
    if (d_ > 2) {
      gradPhysicalZWeighted_ = scalar_view("gradPhysicalZWeighted_", c_, f_, p_);
    }
    divPhysicalWeighted_  = scalar_view("divPhysicalWeighted_", c_, f_, p_);
    gradgradMats_         = scalar_view("gradgradMats_", c_, f_, f_);
    valvalMats_           = scalar_view("valvalMats_", c_, f_, f_);
    cubPointsPhysical_    = scalar_view("cubPointsPhysical_", c_, p_, d_);
    dofPoints_            = scalar_view("dofPoints_", f_, d_);

    /*** START: Fill multidimensional arrays. ***/

    // Compute cubature points and weights.
    cubature_->getCubature(cubPointsSubcell_, cubWeights_);

    // Compute reference basis value and gradient.
    ct::mapToReferenceSubcell(cubPoints_, cubPointsSubcell_, sd_, sideId_, *cellTopo_);
    basis_->getValues(gradReference_, cubPoints_, Intrepid2::OPERATOR_GRAD);       // evaluate grad operator at cubature points
    basis_->getValues(valReference_, cubPoints_, Intrepid2::OPERATOR_VALUE);       // evaluate value operator at cubature points

    // Compute cell Jacobian matrices, its inverses and determinants.
    ct::setJacobian(cellJac_, cubPoints_, cellNodes_, *cellTopo_);  // compute cell Jacobians
    ct::setJacobianInv(cellJacInv_, cellJac_);                      // compute inverses of cell Jacobians
    ct::setJacobianDet(cellJacDet_, cellJac_);                      // compute determinants of cell Jacobians

    // Compute weighted cell measure.
    scalar_view scratch("scratch", c_*p_*d_*d_);
    if (d_ == 2) {
      fst::computeEdgeMeasure(cellWeightedMeasure_, cellJac_, cubWeights_, sideId_, *cellTopo_, scratch);
    }
    else if (d_ == 3) {
      fst::computeFaceMeasure(cellWeightedMeasure_, cellJac_, cubWeights_, sideId_, *cellTopo_, scratch);
    }

    // Transform reference values into physical space.
    fst::HGRADtransformVALUE(valPhysical_,valReference_);

    // Multiply with weighted measure to get weighted values in physical space.
    fst::multiplyMeasure(valPhysicalWeighted_, cellWeightedMeasure_, valPhysical_);

    // Transform reference gradients into physical space.
    fst::HGRADtransformGRAD(gradPhysical_, cellJacInv_, gradReference_);

    // Multiply with weighted measure to get weighted gradients in physical space.
    fst::multiplyMeasure(gradPhysicalWeighted_, cellWeightedMeasure_, gradPhysical_);

    // Extract individual (x, y, z) components of the weighted gradients in physical space.
    for (int c=0; c<c_; ++c) {
      for (int f=0; f<f_; ++f) {
        for (int p=0; p<p_; ++p) {
          gradPhysicalX_(c,f,p) = gradPhysical_(c,f,p,0);
          gradPhysicalXWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,0);
          if (d_ > 1) {
          gradPhysicalY_(c,f,p) = gradPhysical_(c,f,p,1);
          gradPhysicalYWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,1);
          }
          if (d_ > 2) {
          gradPhysicalZ_(c,f,p) = gradPhysical_(c,f,p,2);
          gradPhysicalZWeighted_(c,f,p) = gradPhysicalWeighted_(c,f,p,2);
          }
        }
      }
    }

    // Build divergence in physical space.
    rst::add(divPhysical_, gradPhysicalX_);
    if (d_ > 1) {
      rst::add(divPhysical_, gradPhysicalY_);
    }
    if (d_ > 2) {
      rst::add(divPhysical_, gradPhysicalZ_);
    }

    // Multiply with weighted measure to get weighted divegence in physical space.
    fst::multiplyMeasure(divPhysicalWeighted_, cellWeightedMeasure_, divPhysical_);

    // Compute stiffness matrices.
    fst::integrate(gradgradMats_, gradPhysical_, gradPhysicalWeighted_);

    // Compute mass matrices.
    fst::integrate(valvalMats_, valPhysical_, valPhysicalWeighted_);

    // Map reference cubature points to cells in physical space.
    ct::mapToPhysicalFrame(cubPointsPhysical_,
                                                  cubPoints_,
                                                  cellNodes_,
                                                  *cellTopo_);

    // Compute local degrees of freedom on reference cell sides.
    int numSides = cellTopo_->getSideCount();
    for (int i=0; i<numSides; ++i) {
      sideDofs_.push_back(computeBoundaryDofs(i));
    }

    // Get coordinates of DOFs in reference cell.
      basis->getDofCoords(dofPoints_);

    /*** END: Fill multidimensional arrays. ***/

  }

  /** \brief  Returns cell Jacobian matrices at cubature points.
  */
  scalar_view J() const {
    return cellJac_;
  }

  /** \brief  Returns inverses of cell Jacobians at cubature points.
  */
  scalar_view invJ() const {
    return cellJacInv_;
  }

  /** \brief  Returns determinants of cell Jacobians at cubature points.
  */
  scalar_view detJ() const {
    return cellJacDet_;
  }

  /** \brief  Returns values of FE basis at cubature points in reference space.
  */
  scalar_view Nref() const {
    return valReference_;
  }

  /** \brief  Returns gradients of FE basis at cubature points in reference space.
  */
  scalar_view gradNref() const {
    return gradReference_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space.
  */
  scalar_view N() const {
    return valPhysical_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  scalar_view NdetJ() const {
    return valPhysicalWeighted_;
  }

  /** \brief  Returns gradient of FE basis at cubature points in physical space.
  */
  scalar_view gradN() const {
    return gradPhysical_;
  }

  /** \brief  Returns gradient of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  scalar_view gradNdetJ() const {
    return gradPhysicalWeighted_;
  }

  /** \brief  Returns divergence of FE basis at cubature points in physical space.
  */
  scalar_view divN() const {
    return divPhysical_;
  }

  /** \brief  Returns divergence of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  scalar_view divNdetJ() const {
    return divPhysicalWeighted_;
  }

  /** \brief  Returns x, y or z component of the gradient of FE basis at
              cubature points in physical space.

      \param  coord    [in]   - coordinate index (x=0, y=1, z=2)
  */
  scalar_view DND(const int & coord) const {
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
  scalar_view DNDdetJ(const int & coord) const {
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
  scalar_view stiffMat() const {
    return gradgradMats_;
  }

  /** \brief  Returns mass matrices on cells.
  */
  scalar_view massMat() const {
    return valvalMats_;
  }

  /** \brief Returns cubature points on cells in physical space.
  */
  scalar_view cubPts() const {
    return cubPointsPhysical_;
  }

  /** \brief Builds FE value interpolant and evaluates it at cubature
             points in physical space.
  */
  void evaluateValue(scalar_view fVals,
                     const scalar_view inCoeffs) const {
    fst::evaluate(fVals, inCoeffs, valPhysical_);
  }

  /** \brief Builds FE gradient interpolant and evaluates it at cubature
             points in physical space.
  */
  void evaluateGradient(scalar_view fGrads,
                        const scalar_view inCoeffs) const {
    fst::evaluate(fGrads, inCoeffs, gradPhysical_);
  }

  /** \brief Computes integral of the product or dot-product of interpolated
             FE fields f1 and f2, indexed by (C,P), for values, or (C,P,D),
             for gradients.
  */
  void computeIntegral(scalar_view integral,
                       const scalar_view f1,
                       const scalar_view f2,
                       const bool sumInto = false) const {
    scalar_view f2Weighted;
    if(f2.rank() == 3) {
      f2Weighted = scalar_view("f2Weighted", f2.extent(0), f2.extent(1), f2.extent(2));
    } else { //rank = 2
      f2Weighted = scalar_view("f2Weighted", f2.extent(0), f2.extent(1));
    }
    fst::scalarMultiplyDataData(f2Weighted, cellWeightedMeasure_, f2);  // multiply with weighted measure
                                                              
    fst::integrate(integral, f1, f2Weighted, sumInto);                  // compute integral of f1*f2
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

    auto tagToId = basis_->getAllDofOrdinal();
    std::vector<std::vector<std::vector<int> > > tagToIdCompact;
    int scdim = tagToId.extent_int(0);
    tagToIdCompact.resize(scdim);
    for (int i=0; i<scdim; ++i) {
      for (int j=0; j<tagToId.extent_int(1); ++j) {
        std::vector<int> ids;
        for (int k=0; k<tagToId.extent_int(2); ++k) {
          //std::cout << "\n  i=" << i << "  j=" << j << "  k=" << k << "  id=" << tagToId[i][j][k];
          if (tagToId(i,j,k) != -1) {
            ids.push_back(tagToId(i,j,k));
          }
        }
        tagToIdCompact[i].push_back(ids);
        //int dofordcompact = tagToIdCompact[i][j].size();
        //for (int k=0; k<dofordcompact; ++k) {
        //  std::cout << "\n  i=" << i << "  j=" << j << "  k=" << k << "  id=" << tagToIdCompact[i][j][k];
        //}
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
  void computeDofCoords(scalar_view dofCoords,
                        const scalar_view cellNodes) const {
    // Map reference DOF locations to physical space.
    ct::mapToPhysicalFrame(dofCoords, dofPoints_, cellNodes, *cellTopo_);
  }

}; // FE

#endif
