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

template<class Scalar>
LocalForm0<Scalar>::LocalForm0(Teuchos::RCP<Basis<Scalar> > basis,
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



template<class Scalar>
const VarContainer<Scalar> & LocalForm0<Scalar>::getOperator(const EOperator  primOp,
                                                             const int        subDim,
                                                             const int        subCellId) {

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

  if ((int)basisVals_[primOp].size() == 0) {
    basisVals_[primOp].resize(cubature_.size());
  }
  if ((int)basisVals_[primOp][dimIndex].size() == 0) {
    ((basisVals_[primOp])[dimIndex]).resize(cubature_[dimIndex].size());
  }
  if ((int)basisVals_[primOp][dimIndex][subCellId].getSize() == 0) {
    basisVals_[primOp][dimIndex][subCellId].reset(numCubPoints_[dimIndex][subCellId],
                                                  basisNumDofs_,
                                                  FIELD_FORM_0,
                                                  primOp,
                                                  myCellDim);
    basis_->getValues(basisVals_[primOp][dimIndex][subCellId], cubPoints_[dimIndex][subCellId], primOp);
  }

  return basisVals_[primOp][dimIndex][subCellId];
}



template<class Scalar>
void LocalForm0<Scalar>::getOperator(VarContainer<Scalar> &                  outputValues,
                                     const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                     const EOperator                         primOp) {
  basis_->getValues(outputValues, inputPoints, primOp);
}



template<class Scalar>
void LocalForm0<Scalar>::getOperator(LexContainer<Scalar> &                  outputValues,
                                     const Teuchos::Array<Point<Scalar> > &  inputPoints,
                                     const EOperator                         primOp,
                                     const Cell<Scalar> &                    cell) {
}



template<class Scalar>
void LocalForm0<Scalar>::getOperator(LexContainer<Scalar> &          outputValues,
                                     const EOperator                 leftOp,
                                     const EOperator                 rightOp,
                                     MultiCell<Scalar> &             mCell,
                                     const Teuchos::Array<Scalar> &  inputData,
                                     const EDataFormat               inputFormat,
                                     const EIntegrationDomain        intDomain) {
}



template<class Scalar>
void LocalForm0<Scalar>::getOperator(LexContainer<Scalar> &      outputValues,
                                     const EOperator             leftOp,
                                     const EOperator             rightOp,
                                     MultiCell <Scalar> &        mCell,
                                     const EIntegrationDomain    intDomain) {
  LexContainer<Scalar> leftValues;
  LexContainer<Scalar> rightValues;
  fillLeft(leftValues, leftOp, rightOp, *this, mCell, intDomain);
  fillRight(rightValues, rightOp, leftOp, *this, mCell, intDomain);
  integrate(outputValues, leftValues, rightValues); 
}



template<class Scalar>
void LocalForm0<Scalar>::getOperator(LexContainer<Scalar> &           outputValues,
                                     const EOperator                  leftOp,
                                     const EOperator                  rightOp,
                                     const LocalField<Scalar> &       rightOpField,
                                     MultiCell<Scalar> &              mCell,
                                     const Teuchos::Array<Scalar> &   inputData,
                                     const EDataFormat                inputFormat,
                                     const EIntegrationDomain         intDomain) {
}



template<class Scalar>
void LocalForm0<Scalar>::fillLeft(LexContainer<Scalar> &           leftValues,
                                  const EOperator                  leftOp,
                                  const EOperator                  rightOp,
                                  const LocalField<Scalar> &       rightOpField,
                                  MultiCell<Scalar> &              mCell,
                                  const EIntegrationDomain         intDomain) {

  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }

  switch(leftOp) {                        // Left operator: OPERATOR_VALUE, OPERATOR_GRAD, etc.

   case OPERATOR_VALUE: {

     switch(rightOp) {                      // Right operator:  OPERATOR_VALUE, OPERATOR_GRAD, etc.

       case OPERATOR_VALUE: {

         switch(rightOpField.getFieldType()) { // Right operator field type:  FIELD_FORM_0, etc.

           case FIELD_FORM_0: {
             if (intDomain == INTEGRATION_DOMAIN_CELL) {
               Teuchos::Array<int> iRange(3);
               iRange[0] = 1000;
               iRange[1] = numCubPoints_[0][0];
               iRange[2] = basisNumDofs_;
               leftValues.resize(iRange);
               VarContainer<Scalar> funcVals = getOperator(leftOp, mCell.getMyCellDim(), 0);
               Teuchos::Array<int> miLV(3);
               Teuchos::Array<int> miFV(2);
               for (int cl=0; cl<iRange[0]; cl++) {
                 miLV[0] = cl;
                 for (int qp=0; qp<iRange[1]; qp++) {
                   Matrix<Scalar> jacMat(mCell.getMyCellDim());
                   jacMat = mCell.jacobian(cl, cubPoints_[0][0][qp]);
                   Scalar jacDet = jacMat.det();
                   miLV[1] = qp; miFV[0] = qp;
                   for (int bf=0; bf<iRange[2]; bf++) {
                     miLV[2] = bf; miFV[1] = bf;
                     leftValues.setValue(jacDet*funcVals.getValue(miFV), miLV);
                   }
                 }
               }
             }
           }
           break; //end case FIELD_FORM_0

           default:
             TEST_FOR_EXCEPTION((rightOpField.getFieldType() != FIELD_FORM_0),
                                std::invalid_argument,
                                ">>> ERROR (LocalForm0): Invalid right field type!");

         } // end switch(rightOpField.getFieldType())

       }
       break; // end case OPERATOR VALUE

       default:
         TEST_FOR_EXCEPTION((rightOp != OPERATOR_VALUE),
                            std::invalid_argument,
                            ">>> ERROR (LocalForm0): Invalid right operator!");

     } // end switch(rightOp)

   }
   break; // end case OPERATOR_VALUE

   case OPERATOR_GRAD: {
   }
   break; // end case OPERATOR_GRAD

   default:
     TEST_FOR_EXCEPTION((leftOp != OPERATOR_VALUE) && (leftOp != OPERATOR_GRAD),
                        std::invalid_argument,
                        ">>> ERROR (LocalForm0): Invalid left operator!");

  } // end switch(leftOp)

} // end fillLeft



template<class Scalar>
void LocalForm0<Scalar>::fillRight(LexContainer<Scalar> &           rightValues,
                                   const EOperator                  rightOp,
                                   const EOperator                  leftOp,
                                   const LocalField<Scalar> &       rightOpField,
                                   MultiCell<Scalar> &              mCell,
                                   const EIntegrationDomain         intDomain) {

  if (mCell.getAtlasStatus() == STATUS_UNDEFINED) {
    mCell.setAtlas();
  }

  switch(rightOp) {                       // Right operator: OPERATOR_VALUE, OPERATOR_GRAD, etc.

    case OPERATOR_VALUE: {

      switch(leftOp) {                       // Left operator:  OPERATOR_VALUE, OPERATOR_GRAD, etc.

        case OPERATOR_VALUE: {

          switch(rightOpField.getFieldType()) { // Right operator field type:  FIELD_FORM_0, etc.

            case FIELD_FORM_0: {
              if (intDomain == INTEGRATION_DOMAIN_CELL) {
                Teuchos::Array<int> iRange(3);
                iRange[0] = 1000;
                iRange[1] = numCubPoints_[0][0];
                iRange[2] = basisNumDofs_;
                rightValues.resize(iRange);
                VarContainer<Scalar> funcVals = getOperator(leftOp, mCell.getMyCellDim(), 0);
                Teuchos::Array<int> miRV(3);
                Teuchos::Array<int> miFV(2);
                for (int cl=0; cl<iRange[0]; cl++) {
                  miRV[0] = cl;
                  for (int qp=0; qp<iRange[1]; qp++) {
                    miRV[1] = qp; miFV[0] = qp;
                    for (int bf=0; bf<iRange[2]; bf++) {
                      miRV[2] = bf; miFV[1] = bf;
                      rightValues.setValue(cubWeights_[0][0][qp]*funcVals.getValue(miFV), miRV);
                    }
                  }
                }
              }
            }
            break; //end case FIELD_FORM_0

            default:
              TEST_FOR_EXCEPTION((rightOpField.getFieldType() != FIELD_FORM_0),
                                 std::invalid_argument,
                                 ">>> ERROR (LocalForm0): Invalid right field type!");

          } // end switch(rightOpField.getFieldType())

        }
        break; // end case OPERATOR VALUE

        default:
          TEST_FOR_EXCEPTION((rightOp != OPERATOR_VALUE),
                             std::invalid_argument,
                             ">>> ERROR (LocalForm0): Invalid right operator!");

      } // end switch(leftOp)

    }
    break; // end case OPERATOR_VALUE

    case OPERATOR_GRAD: {
    }
    break; // end case OPERATOR_GRAD

    default:
      TEST_FOR_EXCEPTION((leftOp != OPERATOR_VALUE) && (leftOp != OPERATOR_GRAD),
                         std::invalid_argument,
                         ">>> ERROR (LocalForm0): Invalid left operator!");

  } // end switch(rightOp)

} // end fillRight



template<class Scalar>
void LocalForm0<Scalar>::integrate(LexContainer<Scalar> &        outputValues,
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
  iRange[0] = 1000;
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
                  tmpVal += leftValues.getValue(miLeft)*rightValues.getValue(miRight);
                }
                break;
                case 4: { // vector fields
                  int iBound3 = leftValues.getIndexBound(3);
                  for (int iVec=0; iVec<iBound3; iVec++) {
                    miLeft[3] = iVec; miRight[3] = iVec;
                    tmpVal += leftValues.getValue(miLeft)*rightValues.getValue(miRight);
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
            outputValues.setValue(tmpVal, miOut);
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




}// end namespace Intrepid
