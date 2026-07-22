// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  fe_curl.hpp
    \brief Given a set of cells with geometric node information, sets up
           data structures used in finite element integration in HCURL.
*/

#ifndef PDEOPT_FE_CURL_HPP
#define PDEOPT_FE_CURL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_Cubature.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"

template <class Real>
class FE_CURL {

private:

  const ROL::Ptr<Intrepid::FieldContainer<Real> > cellNodes_;                            // coordinates of the cell nodes
  ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basis_;              // Intrepid basis
  const ROL::Ptr<Intrepid::Cubature<Real, Intrepid::FieldContainer<Real> > > cubature_;  // Intrepid cubature (quadrature, integration) rule
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
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPoints_;             // points of the cubature rule on the reference cell
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubWeights_;            // weights of the cubature rule on the reference cell
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsSubcell_;      // cubature points on the side reference cell
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubWeightsSubcell_;     // cubature weights on the side reference cell
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJac_;               // cell Jacobian matrices
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJacInv_;            // inverses of cell Jacobians
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellJacDet_;            // determinants of cell Jacobians
  ROL::Ptr<Intrepid::FieldContainer<Real> > cellWeightedMeasure_;   // cell measure (Jacobian determinant) multiplied by the cubature weights
  ROL::Ptr<Intrepid::FieldContainer<Real> > valReference_;          // value of FE basis in reference space
  ROL::Ptr<Intrepid::FieldContainer<Real> > curlReference_;         // curl of FE basis in reference space
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysical_;           // value of FE basis in physical space
  ROL::Ptr<Intrepid::FieldContainer<Real> > curlPhysical_;          // curl of FE basis in physical space
  ROL::Ptr<Intrepid::FieldContainer<Real> > valPhysicalWeighted_;   // value of FE basis in physical space multiplied by weighted cell measure
  ROL::Ptr<Intrepid::FieldContainer<Real> > curlPhysicalWeighted_;  // curl of FE basis in physical space multiplied by weighted cell measure
  ROL::Ptr<Intrepid::FieldContainer<Real> > curlcurlMats_;          // cell curl-curl matrices
  ROL::Ptr<Intrepid::FieldContainer<Real> > valvalMats_;            // cell val-val matrices
  ROL::Ptr<Intrepid::FieldContainer<Real> > cubPointsPhysical_;     // cubature points on the physical cells
  ROL::Ptr<Intrepid::FieldContainer<Real> > dofPoints_;             // degree of freedom points on the reference cell

public:

  FE_CURL(const ROL::Ptr<Intrepid::FieldContainer<Real> >                            & cellNodes,
          const ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > >    & basis,
          const ROL::Ptr<Intrepid::Cubature<Real, Intrepid::FieldContainer<Real> > > & cubature) :
    cellNodes_(cellNodes), basis_(basis), cubature_(cubature), sideId_(-1) {

    // Get base cell topology from basis.
    cellTopo_ = ROL::makePtr<shards::CellTopology>(basis_->getBaseCellTopology());

    // Compute dimensions of multidimensional array members.
    c_  = cellNodes_->dimension(0);
    f_  = basis_->getCardinality();
    p_  = cubature_->getNumPoints();
    d_  = cellTopo_->getDimension();
    sd_ = d_ - 1;

    // Get side subcell topology.
    sideTopo_ = ROL::nullPtr;

    // Allocate multidimensional arrays.
    cubPoints_               = ROL::makePtr<Intrepid::FieldContainer<Real>>(p_, d_);
    cubWeights_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(p_);
    cellJac_                 = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, p_, d_, d_);
    cellJacInv_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, p_, d_, d_);
    cellJacDet_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, p_);
    cellWeightedMeasure_     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, p_);
    valReference_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(f_, p_, d_);
    curlReference_           = ROL::makePtr<Intrepid::FieldContainer<Real>>(f_, p_, d_);
    valPhysical_             = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, p_, d_);
    curlPhysical_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, p_, d_);
    valPhysicalWeighted_     = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, p_, d_);
    curlPhysicalWeighted_    = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, p_, d_);
    curlcurlMats_            = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, f_);
    valvalMats_              = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, f_, f_);
    cubPointsPhysical_       = ROL::makePtr<Intrepid::FieldContainer<Real>>(c_, p_, d_);
    dofPoints_               = ROL::makePtr<Intrepid::FieldContainer<Real>>(f_, d_);

    /*** START: Fill multidimensional arrays. ***/

    // Compute cubature points and weights.
    cubature_->getCubature(*cubPoints_, *cubWeights_);

    // Compute reference basis value and curl.
    basis_->getValues(*curlReference_, *cubPoints_, Intrepid::OPERATOR_CURL);       // evaluate curl operator at cubature points
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
    Intrepid::FunctionSpaceTools::HCURLtransformVALUE<Real>(*valPhysical_,
                                                            *cellJacInv_,
                                                            *valReference_);

    // Multiply with weighted measure to get weighted values in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*valPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *valPhysical_);

    // Transform reference curls into physical space.
    Intrepid::FunctionSpaceTools::HCURLtransformCURL<Real>(*curlPhysical_,
                                                           *cellJac_,
                                                           *cellJacDet_,
                                                           *curlReference_);

    // Multiply with weighted measure to get weighted curls in physical space.
    Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(*curlPhysicalWeighted_,
                                                        *cellWeightedMeasure_,
                                                        *curlPhysical_);

    // Compute curl-curl matrices.
    Intrepid::FunctionSpaceTools::integrate<Real>(*curlcurlMats_,
                                                  *curlPhysical_,
                                                  *curlPhysicalWeighted_,
                                                  Intrepid::COMP_CPP);

    // Compute val-val matrices.
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
    ROL::Ptr<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > > coord_iface =
      ROL::dynamicPtrCast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<Real> > >(basis_);
    if (d_ > 1) {
      coord_iface->getDofCoords(*dofPoints_);
    }

    /*** END: Fill multidimensional arrays. ***/

  }

  /** \brief  Returns cell Jacobian matrices at cubature points.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > J() const {
    return cellJac_;
  }

  /** \brief  Returns inverses of cell Jacobians at cubature points.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > invJ() const {
    return cellJacInv_;
  }

  /** \brief  Returns determinants of cell Jacobians at cubature points.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > detJ() const {
    return cellJacDet_;
  }

  /** \brief  Returns values of FE basis at cubature points in reference space.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > Nref() const {
    return valReference_;
  }

  /** \brief  Returns curls of FE basis at cubature points in reference space.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > curlNref() const {
    return curlReference_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > N() const {
    return valPhysical_;
  }

  /** \brief  Returns value of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > NdetJ() const {
    return valPhysicalWeighted_;
  }

  /** \brief  Returns curl of FE basis at cubature points in physical space.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > curlN() const {
    return curlPhysical_;
  }

  /** \brief  Returns curl of FE basis at cubature points in physical space,
              multiplied by weighted cell measures.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > curlNdetJ() const {
    return curlPhysicalWeighted_;
  }

  /** \brief  Returns curl-curl matrices on cells.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > curlcurlMat() const {
    return curlcurlMats_;
  }

  /** \brief  Returns val-val matrices on cells.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > valvalMat() const {
    return valvalMats_;
  }

  /** \brief Returns cubature points on cells in physical space.
  */
  ROL::Ptr<const Intrepid::FieldContainer<Real> > cubPts() const {
    return cubPointsPhysical_;
  }

  /** \brief Builds FE value interpolant and evaluates it at cubature
             points in physical space.  The input coefficients are
             assumed to be 'signed' (with +/- 1 multiplying the
             coefficient depending on the alignment with the reference
             edge direction).
  */
  void evaluateValue(const ROL::Ptr<Intrepid::FieldContainer<Real> > & fVals,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inCoeffs) const {
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fVals, *inCoeffs, *valPhysical_);
  }

  /** \brief Builds FE value interpolant and evaluates it at cubature
             points in physical space.  The input coefficients are
             combined with edge signs (with +/- 1 multiplying the
             coefficient depending on the alignment with the reference
             edge direction).
  */
  void evaluateValue(const ROL::Ptr<Intrepid::FieldContainer<Real> > & fVals,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inCoeffs,
                     const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inSigns) const {
    Intrepid::FieldContainer<Real> coeffsSigned(*inCoeffs);
    Intrepid::FunctionSpaceTools::applyFieldSigns<Real>(coeffsSigned, *inSigns);
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fVals, coeffsSigned, *valPhysical_);
  }

  /** \brief Builds FE curl interpolant and evaluates it at cubature
             points in physical space.  The input coefficients are
             assumed to be 'signed' (with +/- 1 multiplying the
             coefficient depending on the alignment with the reference
             edge direction).
  */
  void evaluateCurl(const ROL::Ptr<Intrepid::FieldContainer<Real> > & fCurls,
                    const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inCoeffs) const {
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fCurls, *inCoeffs, *curlPhysical_);
  }

  /** \brief Builds FE value interpolant and evaluates it at cubature
             points in physical space.  The input coefficients are
             combined with edge signs (with +/- 1 multiplying the
             coefficient depending on the alignment with the reference
             edge direction).
  */
  void evaluateCurl(const ROL::Ptr<Intrepid::FieldContainer<Real> > & fCurls,
                    const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inCoeffs,
                    const ROL::Ptr<const Intrepid::FieldContainer<Real> > & inSigns) const {
    Intrepid::FieldContainer<Real> coeffsSigned(*inCoeffs);
    Intrepid::FunctionSpaceTools::applyFieldSigns<Real>(coeffsSigned, *inSigns);
    Intrepid::FunctionSpaceTools::evaluate<Real>(*fCurls, coeffsSigned, *curlPhysical_);
  }

  /** \brief Computes integral of the dot product of values or curls of interpolated
             FE fields f1 and f2, indexed by (C,P,D).
  */
  void computeIntegral(const ROL::Ptr<Intrepid::FieldContainer<Real> > & integral,
                       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & f1,
                       const ROL::Ptr<const Intrepid::FieldContainer<Real> > & f2,
                       const bool sumInto = false) const {
    Intrepid::FieldContainer<Real> f2Weighted(*f2);
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real>(f2Weighted,              // multiply with weighted measure
                                                               *cellWeightedMeasure_,
                                                               *f2);
    Intrepid::FunctionSpaceTools::integrate<Real>(*integral,                            // compute integral of f1 . f2
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
  void computeDofCoords(const ROL::Ptr<Intrepid::FieldContainer<Real> > & dofCoords,
                        const ROL::Ptr<const Intrepid::FieldContainer<Real> > & cellNodes) const {
    // Map reference DOF locations to physical space.
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*dofCoords,
                                                  *dofPoints_,
                                                  *cellNodes,
                                                  *cellTopo_);
  }

  void mapRefPointsToPhysical(const ROL::Ptr<Intrepid::FieldContainer<Real> > &px,
                              const ROL::Ptr<const Intrepid::FieldContainer<Real> > &rx) const {
    // Map input reference cell points to cells in physical space.
    Intrepid::CellTools<Real>::mapToPhysicalFrame(*px,
                                                  *rx,
                                                  *cellNodes_,
                                                  *cellTopo_);
  }

}; // FE_CURL

#endif
