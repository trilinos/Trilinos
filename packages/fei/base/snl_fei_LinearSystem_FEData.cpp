/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <math.h>

#include <fei_macros.hpp>

#include <fei_FiniteElementData.hpp>

#include <snl_fei_CommUtils.hpp>
#include <snl_fei_LinearSystem_FEData.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix_Impl.hpp>
#include <snl_fei_Constraint.hpp>
#include <snl_fei_Utils.hpp>

#include <fei_DirichletBCRecord.hpp>
#include <fei_DirichletBCManager.hpp>
#include <fei_BCRecord.hpp>
#include <fei_BCManager.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_FEDataFilter.hpp>

#undef fei_file
#define fei_file "snl_fei_LinearSystem_FEData.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::LinearSystem_FEData::LinearSystem_FEData(fei::SharedPtr<FiniteElementData>& feData,
						  fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
  : commUtilsInt_(),
    localProc_(0),
    numProcs_(1),
    feData_(feData),
    matrixGraph_(matrixGraph),
    dbcManager_(NULL),
    bcManager_(NULL),
    lookup_(NULL)
{
  commUtilsInt_ = matrixGraph->getRowSpace()->getCommUtils();

  localProc_ = commUtilsInt_->localProc();
  numProcs_  = commUtilsInt_->numProcs();
}

//----------------------------------------------------------------------------
snl_fei::LinearSystem_FEData::~LinearSystem_FEData()
{
  delete bcManager_;
  for(unsigned i=0; i<attributeNames_.size(); ++i) {
    delete [] attributeNames_[i];
  }
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::loadEssentialBCs(int numIDs,
			 const int* IDs,
			 int idType,
			 int fieldID,
			 int fieldSize,
			 const double *const *gammaValues,
			 const double *const *alphaValues)
{
  if (bcManager_ == NULL) {
    bcManager_ = new BCManager;
  }

  try {
  bcManager_->addBCRecords(idType, numIDs, IDs, fieldID, fieldSize,
				  gammaValues, alphaValues);
  }
  catch(fei::Exception& exc) {
    FEI_CERR << exc.what()<<FEI_ENDL;
    return(-1);
  }
  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         int offsetIntoField,
                         const double* prescribedValues)
{
  if (bcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager;
  }

  try {
    dbcManager_->addBCRecords(numIDs, idType, fieldID, offsetIntoField,
                              IDs, prescribedValues);
  }
  catch(fei::Exception& exc) {
    FEI_CERR << exc.what()<<FEI_ENDL;
    return(-1);
  }
  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::loadEssentialBCs(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         const int* offsetsIntoField,
                         const double* prescribedValues)
{
  if (bcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager;
  }

  try {
    dbcManager_->addBCRecords(numIDs, idType, fieldID, IDs, offsetsIntoField,
                              prescribedValues);
  }
  catch(fei::Exception& exc) {
    FEI_CERR << exc.what()<<FEI_ENDL;
    return(-1);
  }
  return(0);
}

//----------------------------------------------------------------------------
bool snl_fei::LinearSystem_FEData::eqnIsEssentialBC(int globalEqnIndex) const
{
  FEI_CERR << "LinearSystem_FEData::eqnIsEssentialBC NOT IMPLEMENTED!!"<<FEI_ENDL;
  return(false);
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_FEData::getEssentialBCs(std::vector<int>& bcEqns,
                                              std::vector<double>& bcVals) const
{
  FEI_CERR << "LinearSystem_FEData::getEssentialBC_Eqns NOT IMPLEMENTED!!"<<FEI_ENDL;
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_FEData::getConstrainedEqns(std::vector<int>& crEqns) const
{
  matrixGraph_->getConstrainedIndices(crEqns);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::loadComplete(bool applyBCs,
                                               bool globalAssemble)
{
  if (bcManager_ == NULL) {
    bcManager_ = new BCManager;
  }

  int err = 0;
  if (matrix_.get() != NULL && globalAssemble) {
    err = matrix_->gatherFromOverlap();
    if (err != 0) {
      FEI_CERR << "snl_fei::LinearSystem_FEData::loadComplete, ERROR in matrix."
	   << "gatherFromOverlap(), data may be incorrect."<<FEI_ENDL;
    }
  }

  commUtilsInt_->Barrier();

  if (rhs_.get() != NULL && globalAssemble) {
    err = rhs_->gatherFromOverlap();
    if (err != 0) {
      FEI_CERR << "snl_fei::LinearSystem_FEData::loadComplete, ERROR rhs."
        << "gatherFromOverlap(), data may be incorrect."<<FEI_ENDL;
    }
  }

  commUtilsInt_->Barrier();

  CHK_ERR( implementBCs(applyBCs) );

  return(feData_->loadComplete());
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::setBCValuesOnVector(fei::Vector* vector)
{
  ERReturn(-1);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::implementBCs(bool applyBCs)
{
  feiArray<int> essEqns(0, 512);
  feiArray<double> essAlpha(0, 512), essGamma(0, 512);

  feiArray<int> otherEqns(0, 512);
  feiArray<double> otherAlpha(0, 512), otherBeta(0, 512), otherGamma(0, 512);

  fei::SharedPtr<SSMat> localBCeqns(new SSMat);
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph_->getRowSpace();
  int localsize = vecSpace->getNumIndices_Owned();
  fei::Matrix_Impl<SSMat> bcEqns(localBCeqns, matrixGraph_, localsize);

  CHK_ERR( bcManager_->finalizeBCEqns(bcEqns) );

  std::vector<SSMat*>& remote = bcEqns.getRemotelyOwnedMatrix();
  for(unsigned p=0; p<remote.size(); ++p) {
  CHK_ERR( snl_fei::separateBCEqns( *(remote[p]),
				    essEqns, essAlpha, essGamma,
				    otherEqns, otherAlpha, otherBeta, otherGamma) );
  }

  CHK_ERR( bcEqns.gatherFromOverlap() );

  CHK_ERR( snl_fei::separateBCEqns( *(bcEqns.getMatrix()),
				    essEqns, essAlpha, essGamma,
				    otherEqns, otherAlpha, otherBeta, otherGamma) );

  int len = essEqns.length();
  essEqns.resize(len*3);

  int* essEqnsPtr = essEqns.dataPtr();
  double* gammaPtr = essGamma.dataPtr();
  double* alphaPtr = essAlpha.dataPtr();

  for(int i=0; i<len; ++i) {
    int eqn = essEqnsPtr[i];
    essEqnsPtr[i+len] = lookup_->getAssociatedNodeNumber(eqn);
    essEqnsPtr[i+2*len]=lookup_->getOffsetIntoBlkEqn(essEqnsPtr[i+len], eqn);
    double value = gammaPtr[i]/alphaPtr[i];
    gammaPtr[i] = value;
  }

  if (essEqns.length() > 0 && applyBCs) {
    CHK_ERR( feData_->setDirichletBCs(len,
				      essEqnsPtr+len,
				      essEqnsPtr+len*2,
				      gammaPtr) );
  }

  if (otherEqns.length() > 0) {
    ERReturn(-1);
  }

  return(0);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::
loadLagrangeConstraint(int constraintID,
		       const double *weights,
		       double rhsValue)
{
  return(-1);
}

//----------------------------------------------------------------------------
int snl_fei::LinearSystem_FEData::
loadPenaltyConstraint(int constraintID,
		      const double *weights,
		      double penaltyValue,
		      double rhsValue)
{
  return(-1);
}
