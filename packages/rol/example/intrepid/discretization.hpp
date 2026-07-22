// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EXAMPLE_DISCRETIZATION
#define EXAMPLE_DISCRETIZATION

// Sacado Includes
#include "Sacado.hpp"

// Intrepid Includes
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_ArrayTools.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_Utils.hpp"

// Shards Includes
#include "Shards_CellTopology.hpp"


template<class Real> 
class Discretization {

  typedef Intrepid::FunctionSpaceTools                 FST;
  typedef Intrepid::FieldContainer<Real>               FC;
  typedef Intrepid::Cubature<Real>                     Cubature;
  typedef Intrepid::CellTools<Real>                    CellTools;   
  typedef Intrepid::DefaultCubatureFactory<Real>       CubatureFactory;
  typedef Intrepid::Basis_HGRAD_LINE_C1_FEM<Real,FC>   LineBasis;

  template <typename T> using ROL::Ptr = ROL::Ptr<T>;

  private:

    int     numCells_;              // Number of elements
    int     cubDegree_;             // Order of integration cubature
    Real    domainLength_;          // Length of the computational domain     
    int     spaceDim_;              // Number of spatial dimensions    
    int     numNodes_;              // Number of vertices 
    int     numCubPts_;             // Number of cubature points
    int     numFields_;             // Number of fields

    shards::CellTopology cellType_; // Cell Topology
    
    ROL::Ptr<FC> cubPts_;                // Reference cell cubature points
    ROL::Ptr<FC> cubWts_;                // Reference cell cubature weights
    ROL::Ptr<FC> cellNodes_;             // Physical cell vertices
    ROL::Ptr<FC> jacobian_;              // Jacobians of physical cells
    ROL::Ptr<FC> jacInv_;                // Inverses of Jacobians
    ROL::Ptr<FC> jacDet_;                // Determinants of Jacobians
    ROL::Ptr<FC> wtdMeasure_;            // Weighted measure
    ROL::Ptr<FC> valsCubPts_;            // Values at cubature points on reference cell
    ROL::Ptr<FC> gradCubPts_;            // Gradient at cubature points on reference cell
    ROL::Ptr<FC> tranValsCubPts_;        // Transformed values at cubature points 
    ROL::Ptr<FC> tranGradCubPts_;        // Transformed gradient at cubature points
    ROL::Ptr<FC> wtdTranValsCubPts_;     // Weighted transformed values at cubature points   
    ROL::Ptr<FC> wtdTranGradCubPts_;     // Weighted transformed gradient at cubature points 
    ROL::Ptr<FC> physCubPts_;            // Locations of cubature points in physical cells
    ROL::Ptr<FC> massMatrices_;          // Mass matrices for each cells

  public:    

    Discretization(int numCells, int cubDegree, Real domainLength) :
      numCells_(numCells), 
      cubDegree_(cubDegree), 
      domainLength_(domainLength) {

           
 
      // Set the cell topology (intervals) 
      cellType_ = shards::getCellTopologyData<shards::Line<2> >();

      spaceDim_ = cellType_.getDimension();
      numNodes_ = cellType_.getNodeCount();

      CubatureFactory cubFactory; 
 
      // Get the reference cell cubature       
      ROL::Ptr<Cubature> cellCub = cubFactory.create(cellType_,cubDegree_);

      // Get number of cubature points
      numCubPts_ = cellCub->getNumPoints();

      // Create basis
      LineBasis lineBasis;
      numFields_ = lineBasis.getCardinality();

      // Instantiate field containers
      cubPts_            = ROL::makePtr<FC>(numCubPts_,spaceDim_);
      cubWts_            = ROL::makePtr<FC>(numCubPts_);        
      cellNodes_         = ROL::makePtr<FC>(numCells_,numNodes_,spaceDim_);
      jacobian_          = ROL::makePtr<FC>(numCells_,numCubPts_,spaceDim_,spaceDim_);
      jacInv_            = ROL::makePtr<FC>(numCells_,numCubPts_,spaceDim_,spaceDim_);
      jacDet_            = ROL::makePtr<FC>(numCells_,numCubPts_);
      wtdMeasure_        = ROL::makePtr<FC>(numCells_,numCubPts_);
      valsCubPts_        = ROL::makePtr<FC>(numFields_,numCubPts_);
      gradCubPts_        = ROL::makePtr<FC>(numFields_,numCubPts_,spaceDim_);
      tranValsCubPts_    = ROL::makePtr<FC>(numCells_, numFields_, numCubPts_);
      tranGradCubPts_    = ROL::makePtr<FC>(numCells_, numFields_, numCubPts_, spaceDim_);
      wtdTranValsCubPts_ = ROL::makePtr<FC>(numCells_, numFields_, numCubPts_);
      wtdTranGradCubPts_ = ROL::makePtr<FC>(numCells_, numFields_, numCubPts_, spaceDim_);
      physCubPts_        = ROL::makePtr<FC>(numCells_, numCubPts_, spaceDim_);
      massMatrices_      = ROL::makePtr<FC>(numCells_, numFields_, numFields_);

      Real cellLength = domainLength_/numCells_;

      // Set cell vertices (elemental end points)
      for(int i=0; i<numCells_; ++i) {
        (*cellNodes_)(i,0,0) = i*cellLength;
        (*cellNodes_)(i,1,0) = (i+1)*cellLength;
      } 

      // Get the reference cell cubature nodes and weights
      cellCub->getCubature(*cubPts_,*cubWts_);
 
      // Cubature points throughout the domain
      CellTools::mapToPhysicalFrame(*physCubPts_,*cubPts_,*cellNodes_,cellType_);

      // Evaluate the Jacobian, its inverse, and determinant
      CellTools::setJacobian(*jacobian_,*cubPts_,*cellNodes_,cellType_);
      CellTools::setJacobianInv(*jacInv_,*jacobian_);
      CellTools::setJacobianDet(*jacDet_,*jacobian_);
  
      // Evaluate basis functions and their gradients at cubature points
      lineBasis.getValues(*valsCubPts_,*cubPts_,Intrepid::OPERATOR_VALUE);
      lineBasis.getValues(*gradCubPts_,*cubPts_,Intrepid::OPERATOR_GRAD);

      // Get measures and transformed values
      FST::computeCellMeasure<Real>(*wtdMeasure_,*jacDet_,*cubWts_);
      FST::HGRADtransformVALUE<Real>(*tranValsCubPts_,*valsCubPts_);
      FST::HGRADtransformGRAD<Real>(*tranGradCubPts_,*jacInv_,*gradCubPts_); 
      FST::multiplyMeasure<Real>(*wtdTranValsCubPts_,*wtdMeasure_,*tranValsCubPts_);
      FST::multiplyMeasure<Real>(*wtdTranGradCubPts_,*wtdMeasure_,*tranGradCubPts_);
      
      // Store mass matrices
      FST::integrate<Real>(*massMatrices_,*tranValsCubPts_,*wtdTranValsCubPts_,Intrepid::COMP_CPP); 
 
    }        

    
    int getNumCells() {
      return numCells_;
    }

    // Return number of basis functions per element
    int getNumFields() {
      return numFields_; 
    }

    int getSpaceDim() {
      return spaceDim_;
    }
    
    // Return number of cubature points
    int getNumCubPts() {
      return numCubPts_;
    }
    
    ROL::Ptr<FC> getPhysCubPts() {
      return physCubPts_;
    }

    ROL::Ptr<FC> getTransformedVals() {
      return tranValsCubPts_;  
    } 
     
    ROL::Ptr<FC> getTransformedGrad() {
      return tranGradCubPts_;
    }

    ROL::Ptr<FC> getWeightedTransformedVals() {
      return wtdTranValsCubPts_;  
    } 
     
    ROL::Ptr<FC> getWeightedTransformedGrad() {
      return wtdTranGradCubPts_;
    }

    ROL::Ptr<FC> getMassMatrices() {
      return massMatrices_;
    }

};



#endif
