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
const VarContainer<Scalar> & LocalForm0<Scalar,ArrayType>::getOperator(const EOperator  primOp,
                                                                       const int        subDim,
                                                                       const int        subCellId) {

  int myCellDim = MultiCell<Scalar>::getCellDim(basisCell_);
  
  // The first index in basisVals_ is the cell dimension minus the subcell dimension
  int dimIndex  = myCellDim - subDim;
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION(((dimIndex < 0) || (subDim < 0) || (dimIndex >= (int)cubature_.size())),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Invalid subcell dimension or no match for subcell dimension within provided cubature rules!"); 
  TEST_FOR_EXCEPTION(((subCellId < 0) || (subCellId >= (int)cubature_[dimIndex].size())),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Invalid subcell id or no match for subcell id within provided cubature rules!");
#endif
  
  // basisVals_[primOp] is 2-dimensional array of VarContainers whose leading dimension must equal
  // cubature_.size() -  the number of different subcell dimensions for which cubatures are defined
  if ((int)basisVals_[primOp].size() == 0) {
    basisVals_[primOp].resize(cubature_.size());
  }
  
  // basisVals_[primOp][dimIndex] is one-dimensional array of VarContainers whose length must equal
  //  cubature_[dimIndex].size() -  the number of subcells of the specified dimension.
  if ((int)basisVals_[primOp][dimIndex].size() == 0) {
    ((basisVals_[primOp])[dimIndex]).resize(cubature_[dimIndex].size());
  }
  
  if ((int)basisVals_[primOp][dimIndex][subCellId].getSize() == 0) {
    
    // If the VarContainer at basisVals_[primOp][dimIndex][subCellId] is empty, shape it accordingly
    basisVals_[primOp][dimIndex][subCellId].reset(numCubPoints_[dimIndex][subCellId],
                                                  basisNumDofs_,
                                                  FIELD_FORM_0,
                                                  primOp,
                                                  myCellDim);
    
    // Then call getValues from the native basis clas to fill the container.
    basis_ -> getValues(basisVals_[primOp][dimIndex][subCellId], cubPoints_[dimIndex][subCellId], primOp);
  }

  // Otherwise, we simply return const reference to the appropriate VarContainer:
  return basisVals_[primOp][dimIndex][subCellId];
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::transformBasisVals(LexContainer<Scalar> &           transVals,
                                                      const EOperator                  primOp,
                                                      MultiCell<Scalar> &              mCell,
                                                      const bool                       reuseJacobians,
                                                      const EIntegrationDomain         intDomain) {
  // Initializations
  int myCellDim          = mCell.getMyCellDim();
  int numCells           = mCell.getMyNumCells();
  
  // Primitive operator always acts on the native (to the LocalForm0) basis
  switch(primOp) { 
    
    // Admissible operators for LocalForm0: VALUE, GRAD, CURL, D1,...,D10
    case OPERATOR_VALUE: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        int numCubPts = numCubPoints_[0][0];
         
        // for OPERATOR_VALUE transVals is a rank-3 LexContainer where:
        // - 1st index is the number of integration domains, in this case = number of cells in mCell;
        // - 2nd index is the number of cubature points per integration domain;
        // - 3rd index is the number of basis functions in the native basis:
        Teuchos::Array<int> iRange(3);
        iRange[0] = numCells;
        iRange[1] = numCubPts;
        iRange[2] = basisNumDofs_;
        transVals.resize(iRange);
        
        // The appropriate basis values are computed on-demand by this private member function
        VarContainer<Scalar> basisVals = getOperator(primOp, myCellDim, 0);
        
        // Multi-index to access transVals which is rank 3 for value: [cellId][cubPt][bfId]
        Teuchos::Array<int> miTV(3);
        
        // Multi-index to access basisVals which is rank 2 for VALUE: [cubPt][bfId]
        Teuchos::Array<int> miBV(2);
        
        // Loop over cells
        for (int cl = 0; cl < numCells; cl++) {
          miTV[0] = cl;
          
          // Loop over cubature points
          for (int qp = 0; qp < numCubPts; qp++) {
            miTV[1] = qp; 
            miBV[0] = qp;
            
            // Loop over basis functions
            for (int bf = 0; bf < basisNumDofs_; bf++) {
              miTV[2] = bf; 
              miBV[1] = bf;
              
              //transVals.setValue(basisVals.getValue(miBV), miTV);
              //Use this if performance is critical: 
              transVals.getData()[(cl*iRange[1]+qp)*iRange[2]+bf] =  basisVals.getData()[qp*iRange[2]+bf];
            }
          }
        }
      } //if integration domain
    }
      break; // end case OPERATOR VALUE
      
    case OPERATOR_GRAD: 
    case OPERATOR_D1: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        int numCubPts = numCubPoints_[0][0];
        
        // for OPERATOR_GRAD/D1 transVals is a rank-4 LexContainer where:
        // - 1st index < number of integration domains. Range: number of cells in mCell;
        // - 2nd index < number of cubature points per integration domain. Upper bound numCubPoints_[0][0]
        // - 3rd index < number of basis functions in the native basis;
        // - 4th index < dimension of the generating cell type
        Teuchos::Array<int> iRange(4);
        iRange[0] = numCells;
        iRange[1] = numCubPts;
        iRange[2] = basisNumDofs_;
        iRange[3] = myCellDim;
        transVals.resize(iRange);
        
        // The appropriate basis values are computed on-demand by this private member function
        VarContainer<Scalar> basisVals = getOperator(primOp, myCellDim, 0);
        
        // Multi-index to access transVals which is rank 4 for GRAD/D1: [cellId][cubPt][bfId][coord]
        Teuchos::Array<int> miTV(4);
        
        // Multi-index to access basisVals which is rank 3 for GRAD/D1: [cubPt][bfId][coord]
        Teuchos::Array<int> miBV(3);
        
        // This matrix is used if reuseJacobians = false and all Jacobian values are computed on the fly
        Matrix<Scalar> tempJacMat(myCellDim);
        
        // We'll get the type of mapping because affine is handled differently from non-affine
        EMapping cellMapping = MAPPING_MAX;

        // Loop over cells
        for (int cl = 0; cl < numCells; cl++) {
          miTV[0] = cl;
          
          // Get the type the mapping because affine is handled differently from non-affine
          cellMapping = mCell.getCellMappingType(cl);

          // Loop over cubature points
          for (int qp = 0; qp < numCubPts; qp++) {
            miTV[1] = qp; 
            miBV[0] = qp;
            
            // Default is to compute jacobians on the fly
            if( reuseJacobians ){
              // For affine mappings only the first value of the inv. transp. jacobian is stored
              if(cellMapping == MAPPING_AFFINE) {
                // Access only once!
                if( qp == 0 ) {
                  tempJacMat =  mCell.getJacobianTInv(cl, myCellDim, 0)[0];
                }
              }
              // For non-affine mappings all values of the inv. transpose jacobian are stored
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
            
            // Loop over basis functions
            for (int bf = 0; bf < basisNumDofs_; bf++) {
              miTV[2] = bf; 
              miBV[1] = bf;
              
              // multi-index of 1st component of the transformed GRAD of basis funct. bf at cub. pt. qp:
              miTV[3] = 0;
              int indexTV = transVals.getEnumeration(miTV);
              
              // multi-index of 1st component of reference GRAD of basis function bf at cub. pt. qp:
              miBV[2] = 0;
              int indexBV = basisVals.getEnumeration(miBV);
              //tempJacMat.multiplyLeft(transVals, indexTV, basisVals, indexBV,iRange[3]);
              //Use this if performance is critical: 
              tempJacMat.multiplyLeft(&transVals.getData()[indexTV], &basisVals.getData()[indexBV]);
            }
          }
        }
      }//if integration domain
    }
      break; // end case OPERATOR_GRAD
      
    default:
      TEST_FOR_EXCEPTION((primOp != OPERATOR_VALUE) && 
                         (primOp != OPERATOR_GRAD),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Invalid primitive operator!");
  } // end switch(primOp)
} // end fillLeft



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::transformBasisVals(LexContainer<Scalar> &           transValues,
                                                      const EOperator                  primOp,
                                                      const LocalField<Scalar> &       primOpField,
                                                      MultiCell<Scalar> &              mCell,
                                                      const bool                       reuseJacobians,
                                                      const EIntegrationDomain         intDomain)
{
  
  // This method acts on an auxiliary LocalField: we will use transformBasisValues() from that field!
  EField primFieldType = primOpField.getFieldType();
  
  // primOpField is of the abstract type LocalField which does not have transformBasisvals method
  // Downcast to the concrete type of the primOpField argument by using it field type:
  switch( primFieldType ) {
    
    case FIELD_FORM_0:
      (dynamic_cast<LocalForm0<Scalar>& > (const_cast<LocalField<Scalar>& >(primOpField) ) ).\
           transformBasisVals(transValues,primOp,mCell,reuseJacobians, intDomain)  ;
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
void LocalForm0<Scalar,ArrayType>::applyWeightedMeasure(LexContainer<Scalar> &         finalVals,
                                                        const LexContainer<Scalar> &   transVals,
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
  int indexFTV           = 0;
  Scalar weightedMeasure = 0.0;
  int myCellDim          = mCell.getMyCellDim();
  int numCells           = transVals.getIndexBound(0);
  int numCubPts          = transVals.getIndexBound(1);
  int numBf              = transVals.getIndexBound(2);
  
  // Temp matrix needed if reuseJacobians = false and all measure data is computed on the fly
  Matrix<Scalar> tempJacMat(myCellDim);
  Scalar         tempJacDet = 0.0;

  // If the cell has an affine mapping we will compute weighted measures once!
  EMapping cellMapping = MAPPING_MAX;

  switch(primOp) {
    case OPERATOR_VALUE: {
      if(intDomain == INTEGRATION_DOMAIN_CELL) {
         
        // Multi-index for transVals and finalVals which are rank 3 for VALUE: [cellId][cubPt][bfId]
        Teuchos::Array<int> miFTV(3);
        
        // Loop over cells: number of cells is the 1st index of the containers
        for (int cl = 0; cl < numCells; cl++) {
          miFTV[0] = cl;
                    
          cellMapping = mCell.getCellMappingType(cl);

          // Loop over cubature points: number of cub. points is the second index of the containers
          for (int qp = 0; qp < numCubPts; qp++) {
            miFTV[1] = qp; 
            
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
            
            // Loop over basis functions: number of basis functions is the third index of the containers
            for (int bf = 0; bf < numBf; bf++) {
              miFTV[2] = bf; 
              
              // Precompute enumeration to avoid computing it twice by each container:
              indexFTV = transVals.getEnumeration(miFTV);
              
              // Apply the weighted measure (Jacobian at cub. pt. times cub. weight) uisng access by enumeration 
              finalVals.setValue(transVals[indexFTV]*weightedMeasure, indexFTV);
              
            }
          }
        }        
      }// if(intDomain)
    }// end case OPERATOR_VALUE
    break;
      
    case OPERATOR_GRAD: 
    case OPERATOR_D1: {
      if (intDomain == INTEGRATION_DOMAIN_CELL) {
        
        // Multi-index for transVals and finalVals which are rank 3 for VALUE: [cellId][cubPt][bfId][coord]
        Teuchos::Array<int> miFTV(4);
        
        // Loop over cells: number of cells is the 1st index of the containers
        for (int cl = 0; cl < numCells; cl++) {
          miFTV[0] = cl;
                    
          cellMapping = mCell.getCellMappingType(cl);

          // Loop over cubature points: number of cub. points is the second index of the containers
          for (int qp = 0; qp < numCubPts; qp++) {
            miFTV[1] = qp; 
            
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
            
            // Loop over basis functions: number of basis functions is the third index of the containers
            for (int bf = 0; bf < numBf; bf++) {
              miFTV[2] = bf; 
              
              //Unroll the dimension loop: get enumeration of the first GRAD component:
              miFTV[3] = 0;
              indexFTV = transVals.getEnumeration(miFTV);
              
              // The 1st component is always there:
              finalVals.setValue(transVals[indexFTV]*weightedMeasure, indexFTV);

              // Set 2nd and 3rd components depending on the space dimension:
              if(myCellDim == 2) finalVals.setValue(transVals[indexFTV + 1]*weightedMeasure, indexFTV + 1);
              if(myCellDim == 3) finalVals.setValue(transVals[indexFTV + 2]*weightedMeasure, indexFTV + 2);
             }
          }
        }
      } // if(intDomain)
    }
    break; // end case OPERTOR_GRAD & D1
      
    default:
      TEST_FOR_EXCEPTION((primOp != OPERATOR_VALUE) && 
                         (primOp != OPERATOR_GRAD),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Invalid primitive operator!");
  }// switch(primOp)
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::integrate(LexContainer<Scalar> &        outputValues,
                                             const LexContainer<Scalar> &  leftValues,
                                             const LexContainer<Scalar> &  rightValues) const {
  int opRank = leftValues.getRank();
#ifdef HAVE_INTREPID_DEBUG
  TEST_FOR_EXCEPTION((opRank != rightValues.getRank()),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Ranks of leftValues and rightValues do not match!");
#endif
  
  int numCells = leftValues.getIndexBound(0);
#ifdef HAVE_INTREPID_DEBUG
  int numRightCells = rightValues.getIndexBound(0);
  TEST_FOR_EXCEPTION((numCells != numRightCells),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of cells in leftValues and rightValues do not agree!");
#endif
  
  int numQps = leftValues.getIndexBound(1);
#ifdef HAVE_INTREPID_DEBUG
  int numRightQps = rightValues.getIndexBound(1);
  TEST_FOR_EXCEPTION((numQps != numRightQps),
                     std::invalid_argument,
                     ">>> ERROR (LocalForm0): Numbers of integration points in leftValues and rightValues do not agree!");
#endif
  
  int numLeftBfs  = leftValues.getIndexBound(2);
  int numRightBfs = rightValues.getIndexBound(2);
  
#ifdef HAVE_INTREPID_DEBUG
  for (int i=3; i<opRank; i++) {
    TEST_FOR_EXCEPTION((leftValues.getIndexBound(i) != rightValues.getIndexBound(i)),
                       std::invalid_argument,
                       ">>> ERROR (LocalForm0): Field sizes in leftValues and rightValues do not agree!");
  }
#endif
  
  Teuchos::Array<int> iRange(3);
  iRange[0] = numCells;
  iRange[1] = numLeftBfs;
  iRange[2] = numRightBfs;
  outputValues.resize(iRange);
  
  switch(compEngine_) {
    
    case COMP_CPP: {
      Teuchos::Array<int> miOut(3);
      Teuchos::Array<int> miLeft;
      Teuchos::Array<int> miRight;
      leftValues.getIndexRange(miLeft);
      rightValues.getIndexRange(miRight);
      
      for (int cl=0; cl<numCells; cl++) {
        miOut[0] = cl; miLeft[0] = cl; miRight[0] = cl;
        
        for (int lbf=0; lbf<numLeftBfs; lbf++) {
          miOut[1] = lbf; miLeft[2] = lbf;
          
          for (int rbf=0; rbf<numRightBfs; rbf++) {
            miOut[2] = rbf; miRight[2] = rbf;
            Scalar tmpVal(0);
            
            for (int qp=0; qp<numQps; qp++) {
              miLeft[1] = qp; miRight[1] = qp;
              
              switch (opRank) {
                case 3: { // scalar fields
                  //tmpVal += leftValues.getValue(miLeft)*rightValues.getValue(miRight);
                  //Use this if performance is critical: 
                  tmpVal += leftValues.getData()[(cl*numQps+qp)*numLeftBfs+lbf]*
                                rightValues.getData()[(cl*numQps+qp)*numRightBfs+rbf];
                }
                break;
                  
                case 4: { // vector fields
                  int iBound3 = leftValues.getIndexBound(3);
                  
                  for (int iVec=0; iVec<iBound3; iVec++) {
                    miLeft[3] = iVec; miRight[3] = iVec;
                    //tmpVal += leftValues.getValue(miLeft)*rightValues.getValue(miRight);
                    //Use this if performance is critical:
                    tmpVal += leftValues.getData()[((cl*numQps+qp)*numLeftBfs+lbf)*iBound3+iVec]*
                                  rightValues.getData()[((cl*numQps+qp)*numRightBfs+rbf)*iBound3+iVec];
                  }
                }
                break;
                  
                case 5: { // tensor fields
                  int iBound3 = leftValues.getIndexBound(3);
                  int iBound4 = leftValues.getIndexBound(4);
                  for (int iTens1=0; iTens1<iBound3; iTens1++) {
                    miLeft[3] = iTens1; miRight[3] = iTens1;
                    for (int iTens2=0; iTens2<iBound4; iTens2++) {
                      miLeft[4] = iTens2; miRight[4] = iTens2;
                      tmpVal += leftValues.getValue(miLeft)*rightValues.getValue(miRight);
                    }
                  }
                }
                break;

                default:
                  TEST_FOR_EXCEPTION(((opRank != 3) && (opRank != 4) && (opRank != 5)),
                                     std::invalid_argument,
                                     ">>> ERROR (LocalForm0): Invalid data rank. Only scalars, vectors, and tensors are supported!");
              }
            }
            //outputValues.setValue(tmpVal, miOut);
            //Use this if performance is critical:
            outputValues.getData()[cl*numLeftBfs*numRightBfs+lbf*numRightBfs+rbf] = tmpVal;
          }
        }
      }
    }
    break;
      
    case COMP_BLAS: {
      int skipL  = numLeftBfs*numQps;
      int skipR  = numRightBfs*numQps;
      int skipOp = numLeftBfs*numRightBfs;
      for (int cl=0; cl<numCells; cl++) {
        Teuchos::BLAS<int, Scalar> myblas;
        double alpha = 1.0;
        double beta  = 0.0;
        myblas.GEMM(Teuchos::NO_TRANS,
                    Teuchos::TRANS,
                    numLeftBfs,
                    numRightBfs,
                    numQps,
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
      
  } // end switch(compEngine_)
  
} // end integrate

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
void LocalForm0<Scalar,ArrayType>::getOperator(VarContainer<Scalar> &                  outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp) {
  basis_->getValues(outputValues, inputPoints, primOp);
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(VarContainer<Scalar> &                  outputValues,
                                               const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                               const EOperator                         primOp,
                                               const Cell<Scalar> &                    cell) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(LexContainer<Scalar> &          outputValues,
                                               const EOperator                 leftOp,
                                               const EOperator                 rightOp,
                                               MultiCell<Scalar> &             mCell,
                                               const Teuchos::Array<Scalar> &  inputData,
                                               const EDataFormat               inputFormat,
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
  LexContainer<Scalar> leftValues;
  LexContainer<Scalar> rightValues;
  
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
void LocalForm0<Scalar,ArrayType>::getOperator(LexContainer<Scalar> &           outputValues,
                                               const EOperator                  leftOp,
                                               const EOperator                  rightOp,
                                               const LocalField<Scalar> &       rightOpField,
                                               MultiCell<Scalar> &              mCell,
                                               const Teuchos::Array<Scalar> &   inputData,
                                               const EDataFormat                inputFormat,
                                               const bool                       reuseJacobians,
                                               const EIntegrationDomain         intDomain) {
}



template<class Scalar, class ArrayType>
void LocalForm0<Scalar,ArrayType>::getOperator(LexContainer<Scalar> &           outputValues,
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
  LexContainer<Scalar> leftValues;
  LexContainer<Scalar> rightValues;
  
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
