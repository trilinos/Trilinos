// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
            
            // Loop over gradient components
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
} // end fillRight 



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



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::integrate(FieldContainer<Scalar> &        outputValues,
                                             const FieldContainer<Scalar> &  leftValues,
                                             const FieldContainer<Scalar> &  rightValues) const {
  
  // leftValues and rightValues can have ranks 3,4,5 and their ranks must be the same!
  // 1st dim. is number of cells; 2nd dim. is number of points.
  int lrRank   = leftValues.getRank();
  int numCells = leftValues.getDimension(0);
  int numQps   = leftValues.getDimension(1);
  
  // outputValues is rank-3 with multi-index (C, L, R). L and R dims are always the last dims in left/rightValues
  int numLeftBfs  = leftValues.getDimension(lrRank  - 1);
  int numRightBfs = rightValues.getDimension(lrRank - 1);
  outputValues.resize(numCells, numLeftBfs, numRightBfs);
  
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION((lrRank != rightValues.getRank()), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Ranks of leftValues and rightValues do not match!");
  int numRightCells = rightValues.getDimension(0);
  TEST_FOR_EXCEPTION((numCells != numRightCells), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of cells in leftValues and rightValues do not agree!");
  int numRightQps = rightValues.getDimension(1);
  TEST_FOR_EXCEPTION((numQps != numRightQps), std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of integration points in leftValues and rightValues do not agree!");
#endif
  
  switch(compEngine_) {
    case COMP_CPP: {
      
      for (int cl=0; cl<numCells; cl++) {
        for (int lbf=0; lbf<numLeftBfs; lbf++) {
          for (int rbf=0; rbf<numRightBfs; rbf++) {
            Scalar tmpVal(0);
            for (int qp=0; qp<numQps; qp++) {
              
              switch (lrRank) {
                case 3: { // scalar fields: multi-index is (C,P,F)
                  tmpVal += leftValues(cl, qp, lbf)*rightValues(cl, qp, rbf);
                }
                break;
                  
                case 4: { // vector fields: multi-index is (C,P,D,F), loop over (D) subindex
                  int vecDim = leftValues.getDimension(2);
                  for (int iVec = 0; iVec < vecDim; iVec++) {
                    tmpVal += leftValues(cl, qp, iVec, lbf)*rightValues(cl, qp, iVec, rbf);
                  }
                }
                break;
                  
                case 5: { // tensor fields: multi-index is (C,P,D,D,F), loop over (D,D) subindex
                  int tenDim0 = leftValues.getDimension(2);
                  int tenDim1 = leftValues.getDimension(3);
                  for (int iTens1 = 0; iTens1 < tenDim0; iTens1++) {
                    for (int iTens2 =0; iTens2 < tenDim1; iTens2++) {
                      tmpVal += leftValues(cl, qp, iTens1, iTens2, lbf)*rightValues(cl, qp, iTens1, iTens2, rbf);
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
      
      // Data size is defined by the rank of the data containers (left and right must have same rank)
      // rank 3: (C,P,F)      -> dataSize = 1 
      // rank 4: (C,P,D,F)    -> dataSize = D     = leftValues.getDimension(2);
      // rank 5: (C,P,D,D,F)  -> dataSize = D*D   = leftValues.getDimension(2)*leftValues.getDimension(3)
      int dataSize = 1;
      if(lrRank == 4) dataSize = leftValues.getDimension(2);
      if(lrRank == 5) dataSize = leftValues.getDimension(2)*leftValues.getDimension(3);

      // GEMM parameters and their values:
      // TRANSA   NO_TRANS
      // TRANSB   TRANS
      // M        #rows(A)               = numLeftBfs
      // N        #cols(B^T)             = numRightBfs
      // K        #cols(A)               = numData = numQps * dataSize
      // ALPHA    1.0
      // A        left data for cell cl  = leftValues.getData()[cl*skipL]
      // LDA      #rows(A)               = numLeftBfs
      // B        right data for cell cl = rightValues.getData()[cl*skipR]
      // LDB      #rows(B)               = numRightBfs
      // BETA     0.0
      // C        result for cell cl     = outputValues.getData()[cl*skipOp]
      // LDC      #rows(C)               = numLeftBfs
      int numData  = numQps*dataSize;       
      int skipL    = numLeftBfs*numData;        // size of the left data chunk per cell
      int skipR    = numRightBfs*numData;       // size of the right data chunk per cell
      int skipOp   = numLeftBfs*numRightBfs;    // size of the output data chunk per cell
      double alpha = 1.0;                       // these are left unchanged by GEMM 
      double beta  = 0.0;

      for (int cl=0; cl<numCells; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    numLeftBfs,
                    numRightBfs,
                    numData,
                    alpha,
                    &leftValues.getData()[cl*skipL],
                    numLeftBfs,
                    &rightValues.getData()[cl*skipR],
                    numRightBfs,
                    beta,
                    &outputValues.getData()[cl*skipOp],
                    numLeftBfs);
      }
    }
    break;
      
    default:
      TEST_FOR_EXCEPTION(((compEngine_ != COMP_CPP) && (compEngine_ != COMP_BLAS)),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Computational engine not defined!");
      
  } // switch(compEngine_)
  
} // integrate

//===========================================================================//
//                                                                           //
//                         Public methods of LocalForm0                      //
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
                          MultiCell<Scalar>::getNumSubcells(basisCell_, MultiCell<Scalar>::getCellDim(basisCell_)-i)),
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



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                             outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp) {
  // here need to copy to the array type
  //basis_->getValues(outputValues, inputPoints, primOp);
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                             outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp,
                                               const Cell<Scalar> &                    cell) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                     outputValues,
                                               const EOperator                 leftOp,
                                               const EOperator                 rightOp,
                                               MultiCell<Scalar> &             mCell,
                                               const ArrayType &               inputData,
                                               const bool                      reuseJacobians,
                                               const EIntegrationDomain        intDomain) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                 outputValues,
                                               const EOperator             leftOp,
                                               const EOperator             rightOp,
                                               MultiCell <Scalar> &        mCell,
                                               const bool                  reuseJacobians,
                                               const EIntegrationDomain    intDomain)
{
  FieldContainer<Scalar> leftValues;
  FieldContainer<Scalar> rightValues;
  
  // If the user had failed to define an atlas we will use the default atlas.
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
  
  // Fill leftValues with the appropriate transformed native basis function values
  transformBasisVals(leftValues,   leftOp,   mCell, reuseJacobians, intDomain);

  // If left and right operators are the same, right values are left values time the measure:
  if(leftOp == rightOp) {
    applyWeightedMeasure(rightValues, leftValues, leftOp, mCell, reuseJacobians, intDomain); 
  }
  
  // Otherwise, rightValues must be computed on their own and then the  measure is applied:
  else {    
    transformBasisVals(rightValues,   rightOp,   mCell, reuseJacobians, intDomain);
    applyWeightedMeasure(rightValues, rightValues, rightOp, mCell, reuseJacobians, intDomain); 
  }
  
  // Dot product of the data assembled in leftValues and rightValues gives the desired integral
  integrate(outputValues, leftValues, rightValues); 
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                      outputValues,
                                               const EOperator                  leftOp,
                                               const EOperator                  rightOp,
                                               const LocalField<Scalar> &       rightOpField,
                                               MultiCell<Scalar> &              mCell,
                                               const ArrayType &                inputData,
                                               const bool                       reuseJacobians,
                                               const EIntegrationDomain         intDomain) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(ArrayType &                      outputValues,
                                               const EOperator                  leftOp,
                                               const EOperator                  rightOp,
                                               const LocalField<Scalar> &       rightOpField,
                                               MultiCell<Scalar> &              mCell,
                                               const bool                       reuseJacobians,
                                               const EIntegrationDomain         intDomain) 
{
#ifdef HAVE_INTREPID_DEBUG
  // The native LocalForm0 and the rightOpFieldLocalField must be instantiated on the same cell type
  TEST_FOR_EXCEPTION( (basisCell_ != rightOpField.getCellType() ), std::invalid_argument,
                      ">>> ERROR (LocalForm0): Right LocalField must be instantiated for the same cell type as the native LocalField!");
#endif
  FieldContainer<Scalar> leftValues;
  FieldContainer<Scalar> rightValues;
  
  // If the user had failed to define an atlas we will use the default atlas.
  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }
  
  // If the option to reuse Jacobian and measure data is selected, precompute and store values
  if(reuseJacobians) {
    if (intDomain == INTEGRATION_DOMAIN_CELL) {
      mCell.initializeMeasures(mCell.getMyCellDim(), 0, cubPoints_[0][0], cubWeights_[0][0]); 
      
      // Check if the native field and the right operator field have cubature sets with matching number
      // of cubature points on the specified integration domain. Note: this test cannot determine
      // whether or not the two cubatures are the same! In this case, dimension of the integration domain
      // is given by mCell.getMyCellDim(). 
#ifdef HAVE_INTREPID_DEBUG
      TEST_FOR_EXCEPTION( (numCubPoints_[0][0] != rightOpField.getNumCubPoints(mCell.getMyCellDim(),0) ),
                          std::invalid_argument,
                          ">>> ERROR (LocalForm0): Right LocalField must be instantiated with the same cubature set as the native LocalField!");
#endif
    }
  }
  
  // Fill leftValues and rightValues with the appropriate transformed native/auxiliary basis function values
  transformBasisVals(   leftValues, leftOp,                  mCell, reuseJacobians, intDomain);
  transformBasisVals(  rightValues, rightOp, rightOpField,   mCell, reuseJacobians, intDomain);
  applyWeightedMeasure(rightValues, rightValues, rightOp, mCell, reuseJacobians, intDomain); 
  
  // Dot product of the data assembled in leftValues and rightValues gives the desired integral
  integrate(outputValues, leftValues, rightValues); 
}



template<class Scalar, class ArrayType>
int    LocalForm0<Scalar,ArrayType>::getNumCubPoints(const int subcellDim,
                                                     const int subcellId) const 
{
#ifdef HAVE_INTREPID_DEBUG
  // Subcell dimension has to be at least 1 (no cubature sets are on nodes) and <= the cell dim
  TEST_FOR_EXCEPTION( !( (0 < subcellDim) && (subcellDim <= MultiCell<Scalar>::getCellDim(basisCell_) ) ),
                      std::invalid_argument,
                      ">>> ERROR (LocalForm_0): Invalid subcell dimension. ");
  
  TEST_FOR_EXCEPTION( !( (0 <= subcellId) && (subcellId < MultiCell<Scalar>::getNumSubcells(basisCell_, subcellDim)) ),
                      std::invalid_argument,
                      ">>> ERROR (MultiCell): Invalid subcell Id. ");
#endif
  
  // The first index is dimIndex = (dimension of the instantiation cell) - (dimension of the subcell)
  // The second index is the subcellId, relative to the instantiation cell template
  return  numCubPoints_[MultiCell<Scalar>::getCellDim(basisCell_) - subcellDim][ subcellId];
}


}// end namespace Intrepid
