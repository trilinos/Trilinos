// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_LocalForm0Def.hpp
    \brief  Definition file for the Intrepid::LocalForm0 class.
    \author Created by P. Bochev and D. Ridzal.
*/


namespace Intrepid {

  //===========================================================================//
  //                                                                           //
  //                         Private methods of LocalForm0                     //
  //                                                                           //
  //===========================================================================//
  
template<class Scalar, class ArrayType>
const FieldContainer<Scalar> & LocalForm0<Scalar,ArrayType>::getOperator(const EOperator  primOp,
                                                                         const int        subDim,
                                                                         const int        subCellId) {

  // Initialize dimIndex = cell dimension - subcell dimension
  int myCellDim = MultiCell<Scalar>::getCellDim(basisCell_);
  int dimIndex  = myCellDim - subDim;
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION(((dimIndex < 0) || (subDim < 0) || (dimIndex >= (int)cubature_.size())),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Invalid subcell dimension or no match for subcell dimension within provided cubature rules!"); 
  TEST_FOR_EXCEPTION(((subCellId < 0) || (subCellId >= (int)cubature_[dimIndex].size())),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Invalid subcell id or no match for subcell id within provided cubature rules!");
#endif

  // basisVals_ is 3-dimensional ragged array of FieldContainer objects with indices [primOp][dimIndex][subcellId]
  // 1st dim. = number of operator types (not all can be applied though) and is set by the constructor.
  // 2nd and 3rd dims. are the same as the 1st and 2nd dims of cubature_ array. This ensures that 
  // basisVals_ has the minimal dimensions needed to store values for all cubatures that were specified in cubature_
   if ((int)basisVals_[primOp].size() == 0) {
    basisVals_[primOp].resize(cubature_.size());
  }
  if ((int)basisVals_[primOp][dimIndex].size() == 0) {
    ((basisVals_[primOp])[dimIndex]).resize(cubature_[dimIndex].size());
  }
  // If the FieldContainer at basisVals_[primOp][dimIndex][subCellId] is empty, resize it accordingly
  if ((int)basisVals_[primOp][dimIndex][subCellId].getSize() == 0) {
    basisVals_[primOp][dimIndex][subCellId].resize(numCubPoints_[dimIndex][subCellId],
                                                   basisNumDofs_,
                                                   FIELD_FORM_0,
                                                   primOp,
                                                   myCellDim);
    
    // Then call getValues from the native basis clas to fill the container.
    basis_ -> getValues(basisVals_[primOp][dimIndex][subCellId], 
                        cubPoints_[dimIndex][subCellId], 
                        primOp);
  }
  // Otherwise, we simply return const reference to the appropriate FieldContainer:
  return basisVals_[primOp][dimIndex][subCellId];
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::transformBasisVals(FieldContainer<Scalar> &         transVals,
                                                      const EOperator                  primOp,
                                                      MultiCell<Scalar> &              mCell,
                                                      const bool                       reuseJacobians,
                                                      const EIntegrationDomain         intDomain) {
  // Initializations
  int myCellDim = mCell.getMyCellDim();
  int numCells  = mCell.getMyNumCells();
  
  // In this method primitive operator always acts on the native (to the LocalForm0) basis
  switch(primOp) { 
    
    // Admissible operators for LocalForm0: VALUE, GRAD, CURL, D1,...,D10
    case OPERATOR_VALUE: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        
        // for OPERATOR_VALUE transVals is a rank-3 FieldContainer with multi-index (C,P,F)
        int numCubPts = numCubPoints_[0][0];
        transVals.resize(numCells, numCubPts, basisNumDofs_);
        
        // The appropriate basis values are computed on-demand by this private member function
        FieldContainer<Scalar> basisVals = getOperator(primOp, myCellDim, 0);
        
        // Loop over cells, cubature points and basis functions:
        for (int cl = 0; cl < numCells; cl++) {
          for (int qp = 0; qp < numCubPts; qp++) {
            for (int bf = 0; bf < basisNumDofs_; bf++) {
              transVals(cl, qp, bf) = basisVals(qp, bf);
            } // F-loop
          } // P-loop
        } // C-loop
      } //if(integration domain)
    } // case OPERATOR_VALUE
      break;
      
    case OPERATOR_GRAD: 
    case OPERATOR_D1: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        
        // for OPERATOR_GRAD and D1 transVals is a rank-4 FieldContainer with multi-index (C,P,D,F)
        int numCubPts = numCubPoints_[0][0];
        transVals.resize(numCells, numCubPts, myCellDim, basisNumDofs_);
        FieldContainer<Scalar> basisVals = getOperator(primOp, myCellDim, 0);
        
        // This matrix is used if reuseJacobians = false and all Jacobian values are computed on the fly
        Matrix<Scalar> tempJacMat(myCellDim);
        
        // We'll get the type of mapping because affine is handled differently from non-affine
        EMapping cellMapping = MAPPING_MAX;

        // Loop over cells, cubature points, partial derivatives and basis functions:
        for (int cl = 0; cl < numCells; cl++) {          
          cellMapping = mCell.getCellMappingType(cl);
          for (int qp = 0; qp < numCubPts; qp++) {
            
            // Default is to compute jacobians on the fly
            if( reuseJacobians ){
              if(cellMapping == MAPPING_AFFINE) {
                // For affine mappings only the first value of the inv. transp. Jacobian is stored: access only once!
                if( qp == 0 ) {
                  tempJacMat =  mCell.getJacobianTInv(cl, myCellDim, 0)[0];
                }
              }
              // For non-affine mappings all values of the inv. transpose Jacobian are stored
              else{
                tempJacMat =  mCell.getJacobianTInv(cl, myCellDim, 0)[qp];
              }
            }
            else{ 
              // For affine mappings compute only the first value of the Jacobian
              if(cellMapping == MAPPING_AFFINE) {
                if( qp == 0 ) {
                  tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                  tempJacMat.transpose();
                  tempJacMat.invert();
                }
              }
              // Otherwise, we have to compute them at every cubature point
              else {
                tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                tempJacMat.transpose();
                tempJacMat.invert();
              }
            }
            
            // Loop over gradient components. transVals dimensions are (C,P,D,F) instead of (C,P,F,D)
            // in order to enable application of BLAS GEMM routine which assumes column-major storage
            for(int dim = 0; dim < myCellDim; dim++) {
              for (int bf = 0; bf < basisNumDofs_; bf++) {
                tempJacMat.rowMultiply(transVals(cl, qp, dim, bf),            // result goes here 
                                       dim,                                   // select row dim of DF^{-T}
                                       basisVals,                             // untransformed basis grad vals
                                       basisVals.getEnumeration(qp, bf, 0));  // begin index of the gradient
              } // F-loop
            } // D-loop
          } // P-loop
        } // C-loop
      } // if(integration domain)
    } // case OPERATOR_GRAD
      break; 
      
    default:
      TEST_FOR_EXCEPTION((primOp != OPERATOR_VALUE) && 
                         (primOp != OPERATOR_GRAD),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Invalid primitive operator!");
  } // switch(primOp)
} // transformBasisVals



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::transformBasisVals(FieldContainer<Scalar> &         transValues,
                                                      const EOperator                  primOp,
                                                      const LocalField<Scalar> &       primOpField,
                                                      MultiCell<Scalar> &              mCell,
                                                      const bool                       reuseJacobians,
                                                      const EIntegrationDomain         intDomain)
{
  // This method acts on an auxiliary LocalField: we will use transformBasisValues() from that field!
  EField primFieldType = primOpField.getFieldType();
  
  // primOpField is of the abstract type LocalField which does not have transformBasisVals method
  // Downcast to the concrete type of the primOpField argument by using its field type:
  switch( primFieldType ) {
    
    case FIELD_FORM_0:
      (dynamic_cast<LocalForm0<Scalar,ArrayType>& > (const_cast<LocalField<Scalar, ArrayType>& >(primOpField) ) ).\
           transformBasisVals(transValues, primOp, mCell, reuseJacobians, intDomain)  ;
      break;
      
    case FIELD_FORM_1:
      break;
      
    case FIELD_FORM_2:
      break;
      
    case FIELD_FORM_3:
      break;
      
    case FIELD_VECTOR:
      break;
      
    case FIELD_TENSOR:
      break;  
      
    default:
      TEST_FOR_EXCEPTION( !( (primFieldType == FIELD_FORM_0) ||
                             (primFieldType == FIELD_FORM_1) ||
                             (primFieldType == FIELD_FORM_2) ||
                             (primFieldType == FIELD_FORM_3) ||
                             (primFieldType == FIELD_VECTOR) ||
                             (primFieldType == FIELD_TENSOR) ),
                          std::invalid_argument,
                          ">>> ERROR (LocalForm0): Invalid auxiliary primitive operator field!");
  }
} // end fillTest 



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::applyWeightedMeasure(FieldContainer<Scalar> &       finalVals,
                                                        const FieldContainer<Scalar> &   transVals,
                                                        const EOperator                primOp,
                                                        MultiCell<Scalar> &            mCell,
                                                        const bool                     reuseJacobians,
                                                        const EIntegrationDomain       intDomain)
{
  // If finalVals is not of the same size as transVals (the input) make it the same size!
  if(finalVals.getSize() != transVals.getSize() ) {
    finalVals.resize(transVals);
  }
  
  // Initializations:
  Scalar weightedMeasure = 0.0;
  int tempEnum           = 0;
  int myCellDim          = mCell.getMyCellDim();
  int numCells           = 0;
  int numCubPts          = 0;
  int numBf              = 0;
  EMapping cellMapping   = MAPPING_MAX;  
  Matrix<Scalar> tempJacMat(myCellDim);
  Scalar      tempJacDet = 0.0;

  switch(primOp) {
    case OPERATOR_VALUE: {
      if(intDomain == INTEGRATION_DOMAIN_CELL) {
        
        // for OPERATOR_VALUE transVals is a rank-3 FieldContainer with multi-index (C,P,F): 
        numCells           = transVals.getDimension(0);     // C-dimension
        numCubPts          = transVals.getDimension(1);     // P-dimension
        numBf              = transVals.getDimension(2);     // F-dimension
        
        for (int cl = 0; cl < numCells; cl++) {
          cellMapping = mCell.getCellMappingType(cl);
          for (int qp = 0; qp < numCubPts; qp++) {
            if(reuseJacobians) {
              weightedMeasure = mCell.getWeightedMeasure(cl, myCellDim, 0)[qp];
            }
            else {
              // If the cell has an affine mapping compute weighted measure once!
              if(cellMapping == MAPPING_AFFINE) {
                if( qp == 0 ) {
                  tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                  tempJacDet = std::abs( tempJacMat.det() );
                }
              }
              else {
                tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                tempJacDet = std::abs( tempJacMat.det() );
              }
              weightedMeasure = tempJacDet*cubWeights_[0][0][qp];
            } // if(! reuseJacobians)
            
            for (int bf = 0; bf < numBf; bf++) {    
              // This saves computing the value enumeration twice
              tempEnum = transVals.getEnumeration(cl, qp, bf);
              finalVals.setValue(transVals[tempEnum]*weightedMeasure,tempEnum);
            } // F-loop
          } // P-loop
        }  // C-loop       
      }// if(intDomain)
    }// case OPERATOR_VALUE
    break;
      
    case OPERATOR_GRAD: 
    case OPERATOR_D1: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        
        // for OPERATOR_GRAD transVals is a rank-4 FieldContainer with multi-index (C,P,D,F): 
        numCells           = transVals.getDimension(0);     // C-dimension
        numCubPts          = transVals.getDimension(1);     // P-dimension
        numBf              = transVals.getDimension(3);     // F-dimension
        
        for (int cl = 0; cl < numCells; cl++) {
          cellMapping = mCell.getCellMappingType(cl);
          for (int qp = 0; qp < numCubPts; qp++) {
            if(reuseJacobians) {
              weightedMeasure = mCell.getWeightedMeasure(cl,myCellDim, 0)[qp];
            }
            else {
              // If the cell has an affine mapping compute weighted measure once!
              if(cellMapping == MAPPING_AFFINE) {
                if( qp == 0 ) {
                  tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                  tempJacDet = std::abs( tempJacMat.det() );
                }
              }
              else {
                tempJacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                tempJacDet = std::abs( tempJacMat.det() );
              }
              weightedMeasure = tempJacDet*cubWeights_[0][0][qp];
            } // if(! reuseJacobians)
            
            for(int dim = 0; dim < myCellDim; dim++){
              for (int bf = 0; bf < numBf; bf++) {
                // This saves computing the value enumeration twice
                tempEnum = transVals.getEnumeration(cl, qp, dim, bf);
                finalVals.setValue(transVals[tempEnum]*weightedMeasure, tempEnum);
              } // F-loop
            } // D-loop
          } // P-loop
        } // C-loop
      } // if(intDomain)
    } // case OPERATOR_GRAD & D1
    break; 
      
    default:
      TEST_FOR_EXCEPTION((primOp != OPERATOR_VALUE) && 
                         (primOp != OPERATOR_GRAD),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Invalid primitive operator!");
  }// switch(primOp)
}

//===========================================================================//
//                                                                           //
//            Public methods of LocalForm0: constructors                     //
//                                                                           //
//===========================================================================//

template<class Scalar, class ArrayType>
LocalForm0<Scalar,ArrayType>::LocalForm0(Teuchos::RCP<Basis<Scalar> > basis,
                                         Teuchos::Array<Teuchos::Array<Teuchos::RCP<Cubature<Scalar> > > > cubature,
                                         ECompEngine compEngine) :
basis_(basis), cubature_(cubature) {
  
  basisCell_        = basis_->getCellType();
  basisCoordSystem_ = basis_->getCoordinateSystem();
  basisType_        = basis_->getBasisType();
  basisDegree_      = basis_->getDegree();
  basisDofTags_     = basis_->getAllLocalDofTags();
  basisNumDofs_     = basisDofTags_.size();
  basisVals_.resize(OPERATOR_MAX);
  
  compEngine_       = compEngine;
  
  
  /***** Prepare the cubPoints_, cubWeights_, and numCubPoints_ arrays. *****/
  
  TEST_FOR_EXCEPTION(( ((int)cubature_.size() != 1) &&
                       ((int)cubature_.size() != MultiCell<Scalar>::getCellDim(basisCell_)) ),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): If face or edge cubatures are used, the number of groups of cubature rules must equal cell dimension!");  
  cubPoints_.resize(cubature_.size());
  
  for (unsigned i=0; i<cubPoints_.size(); i++) {
    if (cubPoints_[i].size() != 0) {
      TEST_FOR_EXCEPTION(((int)cubature_[i].size() !=
                          MultiCell<Scalar>::getCellNumSubcells(basisCell_, MultiCell<Scalar>::getCellDim(basisCell_)-i)),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Number of cubature rules per subcell dimension must equal number of subcells of that dimension!");  
    }
    cubPoints_[i].resize(cubature_[i].size());
  }
  
  cubWeights_.resize(cubPoints_.size());
  for (unsigned i=0; i<cubWeights_.size(); i++) {
    cubWeights_[i].resize(cubPoints_[i].size());
  }
  
  numCubPoints_.resize(cubPoints_.size());
  for (unsigned i=0; i<numCubPoints_.size(); i++) {
    numCubPoints_[i].resize(cubPoints_[i].size());
  }
  
  for (unsigned i=0; i<cubPoints_.size(); i++) {
    for (unsigned j=0; j<cubPoints_[i].size(); j++) {
      // verify that the user has provided required cubature rules
      TEST_FOR_EXCEPTION((is_null(cubature_[i][j])),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Required cubature rule was not initialized!");
      // verify that the cubature cell type matches the subcell type
      TEST_FOR_EXCEPTION((cubature_[i][j]->getCellType() !=
                          MultiCell<Scalar>::getSubcellType(basisCell_, MultiCell<Scalar>::getCellDim(basisCell_)-i, j)),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Cubature cell type does not match subcell type!");
      cubature_[i][j]->getCubature(numCubPoints_[i][j], cubPoints_[i][j], cubWeights_[i][j]);
      // TRANSFORM CUBATURE POINTS ASSOCIATED WITH A LOWER-DIMENSIONAL REFERENCE CELL TO
      // ALL CORRESPONDING SUBCELLS OF THE REFERENCE CELL basisCell_ !!!
    }
  }
}

//===========================================================================//
//                                                                           //
//            Public methods of LocalForm0: getOperator methods              //
//                                                                           //
//===========================================================================//

template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                             outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp) 
{  
  // This method can only be called in the FEM context. Check the basis type of this LocalForm0 object
  if( (BASIS_FEM_DEFAULT <= basisType_) && (basisType_ < BASIS_FVD_DEFAULT) ) {
    
    // Define FieldContainer for getValues method of the basis class to return the values. 
    FieldContainer<Scalar> basisValues;
    basis_ -> getValues(basisValues, inputPoints, primOp);
    
    // basisValues is of rank-2 or rank-3; resize outputValues to match basisvalues dimensions
    if(basisValues.getRank() == 2) {
      outputValues.resize(basisValues.getDimension(0), 
                          basisValues.getDimension(1) );
    }
    else {
      outputValues.resize(basisValues.getDimension(0), 
                          basisValues.getDimension(1), 
                          basisValues.getDimension(2) );
    }
    // For simplicity copy basisValues to outputValues using enumeration instead of nested loops and (i,j,k..) methods.
    int dataSize = basisValues.getSize();
    for(int i = 0; i < dataSize; i++) {
      outputValues[i] = basisValues[i]; 
    }
  }
  // Throw an exception because in the FVD context we need a physical cell in order to compute values.
  else {
    TEST_FOR_EXCEPTION( ( (BASIS_FVD_DEFAULT <= basisType_ ) && (basisType_ <= BASIS_MAX) ), 
                        std::invalid_argument,
                        ">>> ERROR (LocalForm0): This method is undefined for FVD bases.");
  }
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                             outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp,
                                               const Cell<Scalar> &                    cell) 
{
  int numPts = inputPoints.size();
  
  // Reconstruction of values on physical cells depends on the context. These are the FEM contexts:
  if( (BASIS_FEM_DEFAULT <= basisType_) && (basisType_ < BASIS_FVD_DEFAULT) ) {
    
    // In the FEM context we first map physical points to reference points
    Teuchos::Array<Point<Scalar> > refPoints( inputPoints.size() );
    for(int pt = 0; pt < numPts; pt++){
      refPoints[pt] = cell.mapToReferenceCell( inputPoints[pt] ); 
    }    
    // Then we compute reference values using the previous getOperator method. 
    this -> getOperator(outputValues, inputPoints, primOp);
    
    // Finaly, we transform the reference values back to the physical cell.
    switch(primOp) {
      case OPERATOR_VALUE: 
        // Values of LocalForm0 fields don't need any further transformations
        break;
        
        // Gradients of LocalForm0 fields must be properly transformed by DF^{-T}.
      case OPERATOR_GRAD:
      case OPERATOR_D1: {
        
        int myCellDim = cell.getMyCellDim();
        Matrix<Scalar> tempJacMat(myCellDim);
        EMapping cellMapping = cell.getCellMappingType();
        
        for(int pt = 0; pt < numPts; pt++) {
          
          // If the cell is affine, Jacobian can be computed once
          if(cellMapping == MAPPING_AFFINE) {
            if( pt == 0 ) {
              tempJacMat = cell.jacobian(refPoints[pt]);
              tempJacMat.transpose();
              tempJacMat.invert();
            }
          }
          // Otherwise, we have to compute them at every cubature point
          else {
            tempJacMat = cell.jacobian(refPoints[pt]);
            tempJacMat.transpose();
            tempJacMat.invert();
          }
          // Loop over basis functions and transform their gradient components
          for (int bf = 0; bf < basisNumDofs_; bf++) {
              tempJacMat.multiplyLeft(&outputValues(pt, bf, 0),             // result goes here 
                                      &outputValues(pt, bf, 0) );
          } // F-loop
        }//P-loop
      }
        break;
        
      default:
        TEST_FOR_EXCEPTION((primOp != OPERATOR_VALUE) && 
                           (primOp != OPERATOR_GRAD),
                           std::invalid_argument,
                           ">>> ERROR (LocalForm0): Invalid primitive operator!");
    }// switch(primOp)
  }// if(basisType_)
  // And these are the FVD contexts:
  else if( (BASIS_FVD_DEFAULT <= basisType_) && (basisType_ < BASIS_MAX) ){
    TEST_FOR_EXCEPTION( true, std::invalid_argument, ">>> ERROR (LocalForm0): Method not implemented");
  }
  else {
    TEST_FOR_EXCEPTION( !( (BASIS_FEM_DEFAULT <= basisType_) && (basisType_ < BASIS_MAX) ),
                        std::invalid_argument,
                        ">>> ERROR (LocalForm0): This object has an invalid basis type.");
  }
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                     outputValues,
                                               const EOperator                 trialOp,
                                               const EOperator                 testOp,
                                               MultiCell<Scalar> &             mCell,
                                               const ArrayType &               inputData,
                                               const bool                      reuseJacobians,
                                               const EIntegrationDomain        intDomain) 
{
  /*This method computes a discrete operator acting on the native basis set and whose range has the 
  same dimension as the native basis set, i.e., a square matrix with dimensions equal the number
  of the native basis functions. The elements of this matrix are computed by contracting two  
  ("trial" and "test") FieldContainers with properly transformed native basis values, and a third
  container (inputData) with user-supplied data. The admissible combinations of containers are as follows:
  ==================================================================================================
  | trialValues | inputData  | testValues|outputVals. | signature |    discrete operator           |
  =============|============|============|============|===========|=================================
  | (C,P,F)    | (C,P)      | (C,P,F)    |  (C,F,F)   | 323 = 88  | mass + scalar data             |
  | (C,P,D,F)  | (C,P)      | (C,P,D,F)  |  (C,F,F)   | 424 = 114 | stiffness + scalar data        |
  | (C,P,D,F)  | (C,P,D)    | (C,P,D,F)  |  (C,F,F)   | 434 = 119 | stiffness + diag. tensor data  |
  | (C,P,D,F)  | (C,P,D,D)  | (C,P,D,F)  |  (C,F,F)   | 444 = 124 | stiffness + full tensor data   |
  | (C,P,F)    | (C,P,D)    | (C,P,D,F)  |  (C,F,F)   | 334 = 94  | advection + vector data        |
  | (C,P,D,F)  | (C,P,D)    | (C,P,F)    |  (C,F,F)   | 433 = 118 | advection + vector data        |
  ==================================================================================================
  The trial signature is a base (5) value, the test signature is the equivalent decimal value
  C -> number of integration domains
  P -> number of integration points on the reference integration domain
  D -> space dimension, i.e., number of vector and tensor function components
  F -> number of native basis functions
  
  trialOp and testOp must be such that "trial" and "test" FieldContainers and inputData are in
  one of the admissible combinations above. 
  */
  
  // This method is called in the FEM context and requires an atlas. Use default atlas if user failed to define one
  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }
  
  // If the option to reuse Jacobian and measure data is selected, precompute and store values
  if(reuseJacobians) {
    if (intDomain == INTEGRATION_DOMAIN_CELL) {
      mCell.initializeMeasures(mCell.getMyCellDim(), 0, 
                               cubPoints_[0][0], 
                               cubWeights_[0][0]); 
    }
  }
  
  // Containers for properly transformed values of trialOp and testOp applied to the same native basis set:
  FieldContainer<Scalar> trialValues;
  FieldContainer<Scalar> testValues;
  
  // Fill trialValues with the appropriate transformed native basis function values
  transformBasisVals(trialValues,   trialOp,   mCell, reuseJacobians, intDomain);
  
  // If trialOp=testOp, testValues = trialValues times the measure; otherwise compute separately
  if(trialOp == testOp) {
    applyWeightedMeasure(testValues, trialValues, trialOp, mCell, reuseJacobians, intDomain); 
  }
  else {    
    transformBasisVals(testValues,   testOp,   mCell, reuseJacobians, intDomain);
    applyWeightedMeasure(testValues, testValues, testOp, mCell, reuseJacobians, intDomain); 
  }

  // Compute the signature of the three containers
  int signature = trialValues.getRank()*25 + inputData.getRank()*5 + testValues.getRank();
  switch(signature) {
    
    case 88:{        // (scalar)(scalar)(scalar)
      testValues.multiplyScalarData(inputData);
      testValues.contractScalar( outputValues, trialValues, compEngine_);
    }
      break;
      
    case 114: {      // (vector)(scalar)(vector)
      testValues.multiplyScalarData(inputData);
      testValues.contractVector( outputValues, trialValues, compEngine_);
    }
      break;
      
    case 119: {      // (vector)(diag tensor)(vector)
    }
      break;
      
    case 124: {      // (vector)(full tensor)(vector)
    }
      break;
      
    case 94: {       // (scalar)(vector)(vector)
      
    }
      break;
      
    case 118: {      // (vector)(vector)(scalar)
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION( !( (signature == 88 ) ||
                             (signature == 114) ||
                             (signature == 119) ||
                             (signature == 124) ||
                             (signature == 94 ) ||
                             (signature == 118) ), std::invalid_argument,
                          ">>> ERROR (LocalForm0): specified trial and test operators are incompatible with the data array");        
  }// signature
  
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                 outputValues,
                                               const EOperator             trialOp,
                                               const EOperator             testOp,
                                               MultiCell <Scalar> &        mCell,
                                               const bool                  reuseJacobians,
                                               const EIntegrationDomain    intDomain)
{
  /*
   This method computes a discrete operator acting on the native basis set and whose range has the 
   same dimension as the native basis set, i.e., a square NxN matrix where N is the number of the
   native basis functions. The elements of this matrix are computed by contracting the "point" index P
   of a "trial" and a "test" FieldContainer holding properly transformed native basis values. The ranks of 
   the fields in these containers must match and can be either 0 (scalars) or 1 (vectors)
   ================================================================================================
   | field type/rank| trialValues rank/dims |   testValues rank/dims  | outputValues rank/dims    |
   ================================================================================================ 
   |  scalar -> 0   |    3 -> (C,P,N)       |       3 -> (C,P,N)      |       2 -> (C,N,N)        |
   |  vector -> 1   |    4 -> (C,P,D,N)     |       4 -> (C,P,D,N)    |       2 -> (C,N,N)        |
   ================================================================================================ 
        C -> number of integration domains
        P -> number of integration points on the reference integration domain
        N -> number of native basis functions
   
   trialOp and testOp must be such that their action on native basis functions results in fields
   with identical ranks, i.e., either scalars or vectors. 
   */
  
  // This method is called in the FEM context and requires an atlas. Use default atlas if user failed to define one
  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }
  // If the option to reuse Jacobian and measure data is selected, precompute and store values
  if(reuseJacobians) {
    if (intDomain == INTEGRATION_DOMAIN_CELL) {
      mCell.initializeMeasures(mCell.getMyCellDim(), 0, 
                               cubPoints_[0][0], 
                               cubWeights_[0][0]); 
    }
  }
  // Containers for properly transformed values of trialOp and testOp applied to the native basis set:
  FieldContainer<Scalar> trialValues;
  FieldContainer<Scalar> testValues;
  
  // Fill trialValues
  transformBasisVals( trialValues,   trialOp,   mCell, reuseJacobians, intDomain );
  
  //If trialOp==testOp set testValues = trialValues*integration measure, otherwise compute testValues separately 
  if(trialOp == testOp) {
    applyWeightedMeasure( testValues, trialValues, testOp, mCell, reuseJacobians, intDomain ); 
  }
  else {    
    transformBasisVals( testValues,  testOp,   mCell, reuseJacobians, intDomain );
    applyWeightedMeasure( testValues, testValues, testOp, mCell, reuseJacobians, intDomain ); 
  }
  // Contract containers with the native basis values. Ranks will be checked by contract methods.
  if( testValues.getRank() == 3) {
    testValues.contractScalar( outputValues, trialValues, compEngine_ );
  }
  else{
    testValues.contractVector( outputValues, trialValues, compEngine_ );
  }
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                      outputValues,
                                               const EOperator                  trialOp,
                                               const EOperator                  testOp,
                                               const LocalField<Scalar> &       testOpField,
                                               MultiCell<Scalar> &              mCell,
                                               const ArrayType &                inputData,
                                               const bool                       reuseJacobians,
                                               const EIntegrationDomain         intDomain) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                      outputValues,
                                               const EOperator                  trialOp,
                                               const EOperator                  testOp,
                                               const LocalField<Scalar> &       testOpField,
                                               MultiCell<Scalar> &              mCell,
                                               const bool                       reuseJacobians,
                                               const EIntegrationDomain         intDomain) 
{
  /*
   This method computes a discrete operator acting on the native basis set and whose range has the 
   same dimension as the basis set of the specified auxiliary LocalField, i.e., a rectangular 
   AxN matrix where N and A are the numbers of the native and auxiliary basis functions, respectively.
   The elements of this matrix are computed by contracting the "point" index P
   of a "trial" and a "test" FieldContainer holding properly transformed native basis values. The ranks of 
   the fields in these containers must match and can be either 0 (scalars) or 1 (vectors)
   ================================================================================================
   | field type/rank| trialValues rank/dims |   testValues rank/dims  | outputValues rank/dims    |
   ================================================================================================ 
   |  scalar -> 0   |    3 -> (C,P,N)       |       3 -> (C,P,A)      |       2 -> (C,A,N)        |
   |  vector -> 1   |    4 -> (C,P,D,N)     |       4 -> (C,P,D,A)    |       2 -> (C,A,N)        |
   ================================================================================================ 
        C -> number of integration domains
        P -> number of integration points on the reference integration domain
        N -> number of native (trial) basis functions (number of matrix columns)
        A -> number of auxiliary (test) basis functions (number of matrix rows)
   
   trialOp and testOp must be such that their action on native and auxiliary basis functions 
   results in fields with identical ranks, i.e., either scalars or vectors. 
   */
  
#ifdef HAVE_INTREPID_DEBUG
  // The native LocalForm0 and the testOpFieldLocalField must be instantiated on the same cell type
  TEST_FOR_EXCEPTION( (basisCell_ != testOpField.getCellType() ), std::invalid_argument,
                      ">>> ERROR (LocalForm0): Auxiliary (test) LocalField must be instantiated for the same cell type as the native (trial) LocalField!");
#endif
  // This method is called in the FEM context and requires an atlas. Use default atlas if user failed to define one
  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }
  // If the option to reuse Jacobian and measure data is selected, precompute and store values
  if(reuseJacobians) {
    if (intDomain == INTEGRATION_DOMAIN_CELL) {
      mCell.initializeMeasures(mCell.getMyCellDim(), 0, cubPoints_[0][0], cubWeights_[0][0]); 
      
#ifdef HAVE_INTREPID_DEBUG
      /* Cubature sets of the native (trial) and auxiliary (test) fields must be the same for the
      specified integration domain. However, we can only check if they have the same number of points, 
      which does not tell us if they are really the same cubature. For INTEGRATION_DOMAIN_CELL
      dimension of the integration domain is the cell dimension: mCell.getMyCellDim().             
      */
      TEST_FOR_EXCEPTION( (numCubPoints_[0][0] != testOpField.getNumCubPoints(mCell.getMyCellDim(),0) ),
                          std::invalid_argument,
                          ">>> ERROR (LocalForm0): Auxiliary (test) LocalField must be instantiated with the same cubature set as the native LocalField!");
#endif
    }
  }
  // Containers for trial and test values
  FieldContainer<Scalar> trialValues;
  FieldContainer<Scalar> testValues;
  
  // Fill trialValues and testValues with appropriately transformed native/auxiliary basis function values
  transformBasisVals(   trialValues, trialOp,               mCell, reuseJacobians, intDomain);
  transformBasisVals(   testValues,  testOp, testOpField,   mCell, reuseJacobians, intDomain);
  applyWeightedMeasure( testValues,  testValues, testOp,    mCell, reuseJacobians, intDomain); 
  
  // contract methods will check that trial and test values have matching ranks that are either 3 or 4
  if( testValues.getRank() == 3) {
    testValues.contractScalar( outputValues, trialValues, compEngine_);
  }
  else{
    testValues.contractVector( outputValues, trialValues, compEngine_);
  }
}


//===========================================================================//
//                                                                           //
//            Public methods of LocalForm0: getFunctional methods            //
//                                                                           //
//===========================================================================//
template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getFunctional(ArrayType &              outputValues,
                                                 const ArrayType &        inputData,
                                                 const EOperator          testOp,
                                                 MultiCell<Scalar> &      mCell,
                                                 const bool               reuseJacobians,
                                                 const EIntegrationDomain intDomain)
{
  /*
   This method computes a discrete functional acting on the native basis set, i.e., a vector whose 
   length equals the number of native basis functions. The elements of this vector are computed by
   contracting a "trial" trialData container, representing a user-supplied field, with a "test" 
   FieldContainer with the properly transformed native basis values. The ranks of the fields in these 
   containers must match, i.e., they can be either scalar (rank 0) or vector (rank 1) fields. 
   ============================================================================================
   |  field/rank    |  dataArray rank   |     testValues rank      |    outputValues rank     |
   ============================================================================================ 
   |  scalar -> 0   |   2 -> (C,P)      |       3 -> (C,P,F)       |       2 -> (C,F)         |
   |  vector -> 1   |   3 -> (C,P,D)    |       4 -> (C,P,D,F)     |       2 -> (C,F)         |
   ============================================================================================ 
   C -> number of integration domains
   P -> number of integration points on the reference integration domain
   F -> number of native basis functions
   
   testOp must be such that when applied to the native basis the rank of the resulting field matches
   the rank of the user-supplied field in dataArray.
   */
  
  // This method is called in the FEM context and requires an atlas. Use default atlas if user failed to define one
  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }
  // If the option to reuse Jacobian and measure data is selected, precompute and store values
  if(reuseJacobians) {
    if (intDomain == INTEGRATION_DOMAIN_CELL) {
      mCell.initializeMeasures(mCell.getMyCellDim(), 0, 
                               cubPoints_[0][0], 
                               cubWeights_[0][0]); 
    }
  }
  // Container for the transformed native basis values
  FieldContainer<Scalar> testValues;
  
  // Transform the basis values and apply the cubature weights
  transformBasisVals(  testValues,   testOp,   mCell, reuseJacobians, intDomain);
  applyWeightedMeasure(testValues, testValues, testOp, mCell, reuseJacobians, intDomain); 
  
  // Contract native basis values with user-supplied data. Ranks will be checked by contract methods
  if( testValues.getRank() == 3) {
    testValues.contractScalarData( outputValues, inputData, compEngine_);
  }
  else{
    testValues.contractVectorData( outputValues, inputData, compEngine_);
  }
}


//===========================================================================//
//                                                                           //
//            Public methods of LocalForm0: other methods                    //
//                                                                           //
//===========================================================================//

template<class Scalar, class ArrayType>
int    LocalForm0<Scalar,ArrayType>::getNumCubPoints(const int subcellDim,
                                                     const int subcellId) const 
{
#ifdef HAVE_INTREPID_DEBUG
  // Subcell dimension has to be at least 1 (no cubature sets are on nodes) and <= the cell dim
  TEST_FOR_EXCEPTION( !( (0 < subcellDim) && (subcellDim <= MultiCell<Scalar>::getCellDim(basisCell_) ) ),
                      std::invalid_argument,
                      ">>> ERROR (LocalForm_0): Invalid subcell dimension. ");
  
  TEST_FOR_EXCEPTION( !( (0 <= subcellId) && (subcellId < MultiCell<Scalar>::getCellNumSubcells(basisCell_, subcellDim)) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid subcell Id. ");
#endif
  
  // The first index is dimIndex = (dimension of the instantiation cell) - (dimension of the subcell)
  // The second index is the subcellId, relative to the instantiation cell template
  return  numCubPoints_[MultiCell<Scalar>::getCellDim(basisCell_) - subcellDim][ subcellId];
}


// ====== OBSOLETE ===



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::integrate(FieldContainer<Scalar> &        outputValues,
                                             const FieldContainer<Scalar> &  trialValues,
                                             const FieldContainer<Scalar> &  testValues) const {
  
  // This method provides a key functionality needed to compute linear operators and functionals
  // by getOperator and getFunctional methods. These two methods reduce numerical integration to 
  // a dot product of two vectors with "trial" and "test" values, which are stored as multi-dimensional arrays.
  
  // trialValues and testValues can have ranks 3,4,5 and their ranks must be the same!
  // 1st dim. is number of cells; 2nd dim. is number of points.
  int lrRank   = testValues.getRank();
  int numCells = testValues.getDimension(0);
  int numQps   = testValues.getDimension(1);
  
  // outputValues is rank-3 with multi-index (C, L, R). L and R dims are always the last dims in trial/testValues
  int numTrialBfs  = trialValues.getDimension(lrRank  - 1);
  int numTestBfs = testValues.getDimension(lrRank - 1);
  outputValues.resize(numCells, numTrialBfs, numTestBfs);
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION((lrRank != testValues.getRank()), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Ranks of trialValues and testValues do not match!");
  int numTestCells = testValues.getDimension(0);
  TEST_FOR_EXCEPTION((numCells != numTestCells), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of cells in trialValues and testValues do not agree!");
  int numTestQps = testValues.getDimension(1);
  TEST_FOR_EXCEPTION((numQps != numTestQps), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of integration points in trialValues and testValues do not agree!");
#endif
  
  switch(compEngine_) {
    case COMP_CPP: {
      for (int cl=0; cl<numCells; cl++) {
        for (int lbf=0; lbf<numTrialBfs; lbf++) {
          for (int rbf=0; rbf<numTestBfs; rbf++) {
            Scalar tmpVal(0);
            for (int qp=0; qp<numQps; qp++) {
              
              switch (lrRank) {
                case 3: { // scalar fields: multi-index is (C,P,F)
                  tmpVal += trialValues(cl, qp, lbf)*testValues(cl, qp, rbf);
                }
                  break;
                  
                case 4: { // vector fields: multi-index is (C,P,D,F), loop over (D) subindex
                  int vecDim = trialValues.getDimension(2);
                  for (int iVec = 0; iVec < vecDim; iVec++) {
                    tmpVal += trialValues(cl, qp, iVec, lbf)*testValues(cl, qp, iVec, rbf);
                  }
                }
                  break;
                  
                case 5: { // tensor fields: multi-index is (C,P,D,D,F), loop over (D,D) subindex
                  int tenDim0 = trialValues.getDimension(2);
                  int tenDim1 = trialValues.getDimension(3);
                  for (int iTens1 = 0; iTens1 < tenDim0; iTens1++) {
                    for (int iTens2 =0; iTens2 < tenDim1; iTens2++) {
                      tmpVal += trialValues(cl, qp, iTens1, iTens2, lbf)*testValues(cl, qp, iTens1, iTens2, rbf);
                    }
                  }
                }
                  break;
                  
                default:
                  TEST_FOR_EXCEPTION(((lrRank != 3) && (lrRank != 4) && (lrRank != 5)),
                                     std::invalid_argument,
                                     ">>> ERROR (LocalForm0): Invalid data rank. Only scalars, vectors, and tensors are supported!");
              } // switch(lrRank)
            } // P-loop
            
            outputValues(cl, lbf, rbf) = tmpVal;
            
          } // R-loop
        } // L-loop
      } // C-loop      
    }
      break;
      
    case COMP_BLAS: {
      
      // Data size is defined by the rank of the data containers (trial and test must have same rank)
      // rank 3: (C,P,F)      -> dataSize = 1 
      // rank 4: (C,P,D,F)    -> dataSize = D     = trialValues.getDimension(2);
      // rank 5: (C,P,D,D,F)  -> dataSize = D*D   = trialValues.getDimension(2)*trialValues.getDimension(3)
      int dataSize = 1;
      if(lrRank == 4) dataSize = trialValues.getDimension(2);
      if(lrRank == 5) dataSize = trialValues.getDimension(2)*trialValues.getDimension(3);
      
      // GEMM parameters and their values:
      // TRANSA   NO_TRANS
      // TRANSB   TRANS
      // M        #rows(A)               = numTrialBfs
      // N        #cols(B^T)             = numTestBfs
      // K        #cols(A)               = numData = numQps * dataSize
      // ALPHA    1.0
      // A        trial data for cell cl  = trialValues.getData()[cl*skipL]
      // LDA      #rows(A)               = numTrialBfs
      // B        test data for cell cl = testValues.getData()[cl*skipR]
      // LDB      #rows(B)               = numTestBfs
      // BETA     0.0
      // C        result for cell cl     = outputValues.getData()[cl*skipOp]
      // LDC      #rows(C)               = numTrialBfs
      int numData  = numQps*dataSize;       
      int skipL    = numTrialBfs*numData;        // size of the trial data chunk per cell
      int skipR    = numTestBfs*numData;       // size of the test data chunk per cell
      int skipOp   = numTrialBfs*numTestBfs;    // size of the output data chunk per cell
      double alpha = 1.0;                       // these are trial unchanged by GEMM 
      double beta  = 0.0;
      
      for (int cl=0; cl<numCells; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    numTrialBfs,
                    numTestBfs,
                    numData,
                    alpha,
                    &trialValues.getData()[cl*skipL],
                    numTrialBfs,
                    &testValues.getData()[cl*skipR],
                    numTestBfs,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    numTrialBfs);
      }
    }
      break;
      
    default:
      TEST_FOR_EXCEPTION(((compEngine_ != COMP_CPP) && (compEngine_ != COMP_BLAS)),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Computational engine not defined!");
  } // switch(compEngine_)  
                                             } // integrate




}// end namespace Intrepid
