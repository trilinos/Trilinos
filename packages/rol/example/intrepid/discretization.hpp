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

  template <typename T> using RCP = Teuchos::RCP<T>;

  private:

    int     numCells_;              // Number of elements
    int     cubDegree_;             // Order of integration cubature
    Real    domainLength_;          // Length of the computational domain     
    int     spaceDim_;              // Number of spatial dimensions    
    int     numNodes_;              // Number of vertices 
    int     numCubPts_;             // Number of cubature points
    int     numFields_;             // Number of fields

    shards::CellTopology cellType_; // Cell Topology
    
    RCP<FC> cubPts_;                // Reference cell cubature points
    RCP<FC> cubWts_;                // Reference cell cubature weights
    RCP<FC> cellNodes_;             // Physical cell vertices
    RCP<FC> jacobian_;              // Jacobians of physical cells
    RCP<FC> jacInv_;                // Inverses of Jacobians
    RCP<FC> jacDet_;                // Determinants of Jacobians
    RCP<FC> wtdMeasure_;            // Weighted measure
    RCP<FC> valsCubPts_;            // Values at cubature points on reference cell
    RCP<FC> gradCubPts_;            // Gradient at cubature points on reference cell
    RCP<FC> tranValsCubPts_;        // Transformed values at cubature points 
    RCP<FC> tranGradCubPts_;        // Transformed gradient at cubature points
    RCP<FC> wtdTranValsCubPts_;     // Weighted transformed values at cubature points   
    RCP<FC> wtdTranGradCubPts_;     // Weighted transformed gradient at cubature points 
    RCP<FC> physCubPts_;            // Locations of cubature points in physical cells
    RCP<FC> massMatrices_;          // Mass matrices for each cells

  public:    

    Discretization(int numCells, int cubDegree, Real domainLength) :
      numCells_(numCells), 
      cubDegree_(cubDegree), 
      domainLength_(domainLength) {

      using Teuchos::rcp;     
 
      // Set the cell topology (intervals) 
      cellType_ = shards::getCellTopologyData<shards::Line<2> >();

      spaceDim_ = cellType_.getDimension();
      numNodes_ = cellType_.getNodeCount();

      CubatureFactory cubFactory; 
 
      // Get the reference cell cubature       
      RCP<Cubature> cellCub = cubFactory.create(cellType_,cubDegree_);

      // Get number of cubature points
      numCubPts_ = cellCub->getNumPoints();

      // Create basis
      LineBasis lineBasis;
      numFields_ = lineBasis.getCardinality();

      // Instantiate field containers
      cubPts_            = rcp( new FC(numCubPts_,spaceDim_) );
      cubWts_            = rcp( new FC(numCubPts_) );        
      cellNodes_         = rcp( new FC(numCells_,numNodes_,spaceDim_) );
      jacobian_          = rcp( new FC(numCells_,numCubPts_,spaceDim_,spaceDim_) );
      jacInv_            = rcp( new FC(numCells_,numCubPts_,spaceDim_,spaceDim_) );
      jacDet_            = rcp( new FC(numCells_,numCubPts_) );
      wtdMeasure_        = rcp( new FC(numCells_,numCubPts_) );
      valsCubPts_        = rcp( new FC(numFields_,numCubPts_) );
      gradCubPts_        = rcp( new FC(numFields_,numCubPts_,spaceDim_) );
      tranValsCubPts_    = rcp( new FC(numCells_, numFields_, numCubPts_) );
      tranGradCubPts_    = rcp( new FC(numCells_, numFields_, numCubPts_, spaceDim_) );
      wtdTranValsCubPts_ = rcp( new FC(numCells_, numFields_, numCubPts_) );
      wtdTranGradCubPts_ = rcp( new FC(numCells_, numFields_, numCubPts_, spaceDim_) );
      physCubPts_        = rcp( new FC(numCells_, numCubPts_, spaceDim_) );
      massMatrices_      = rcp( new FC(numCells_, numFields_, numFields_) );

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
    
    RCP<FC> getPhysCubPts() {
      return physCubPts_;
    }

    RCP<FC> getTransformedVals() {
      return tranValsCubPts_;  
    } 
     
    RCP<FC> getTransformedGrad() {
      return tranGradCubPts_;
    }

    RCP<FC> getWeightedTransformedVals() {
      return wtdTranValsCubPts_;  
    } 
     
    RCP<FC> getWeightedTransformedGrad() {
      return wtdTranGradCubPts_;
    }

    RCP<FC> getMassMatrices() {
      return massMatrices_;
    }

};



#endif
