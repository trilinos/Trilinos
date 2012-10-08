/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_VectorReducer.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Vector_core.hpp>
#include <fei_Vector.hpp>
#include <fei_CommUtils.hpp>

#undef fei_file
#define fei_file "fei_VectorReducer.cpp"

#include <fei_ErrMacros.hpp>

namespace fei {

//----------------------------------------------------------------------------
VectorReducer::VectorReducer(fei::SharedPtr<fei::Reducer> reducer,
                             fei::SharedPtr<fei::Vector> target,
                             bool isSolutionVector)
  : reducer_(reducer),
    target_(target),
    isSolution_(isSolutionVector)
{
  localProc_ = fei::localProc(target->getVectorSpace()->getCommunicator());
  numProcs_ = fei::numProcs(target->getVectorSpace()->getCommunicator());

  fei::Vector_core* target_core = dynamic_cast<fei::Vector_core*>(target.get());
  if (target_core == NULL) {
    throw std::runtime_error("fei::VectorReducer ERROR, target vector not dynamic_cast-able to fei::Vector_core.");
  }

  fei::SharedPtr<fei::VectorSpace> vecspace = target->getVectorSpace();
  int numEqns = vecspace->getNumIndices_SharedAndOwned();
  std::vector<int> eqns;
  vecspace->getIndices_SharedAndOwned(eqns);

  std::vector<int> overlap;
  for(int i=0; i<numEqns; ++i) {
    if (!reducer->isSlaveEqn(eqns[i])) {
      overlap.push_back(reducer->translateToReducedEqn(eqns[i]));
    }
    else {
      std::vector<int> masters;
      reducer->getSlaveMasterEqns(eqns[i], masters);
      for(unsigned j=0; j<masters.size(); ++j) {
        overlap.push_back(reducer->translateToReducedEqn(masters[j]));
      }
    }
  }

  target_core->setOverlap(overlap.size(), &overlap[0]);
}

//----------------------------------------------------------------------------
VectorReducer::~VectorReducer()
{
}

//----------------------------------------------------------------------------
int VectorReducer::putScalar(double scalar)
{
  return(target_->putScalar(scalar));
}

//----------------------------------------------------------------------------
int VectorReducer::update(double a, const fei::Vector* x, double b)
{
  return(target_->update(a, x, b));
}

//----------------------------------------------------------------------------
int VectorReducer::scatterToOverlap()
{
  return(target_->scatterToOverlap());
}

//----------------------------------------------------------------------------
int VectorReducer::gatherFromOverlap(bool accumulate)
{
  reducer_->assembleReducedVector(isSolution_, *target_);
  target_->setCommSizes();
  return(target_->gatherFromOverlap(accumulate));
}

//----------------------------------------------------------------------------
int VectorReducer::sumIn(int numValues,
                         const int* indices, const double* values,
                         int vectorIndex)
{
  return(reducer_->addVectorValues(numValues, indices, values, true,
                                   isSolution_, vectorIndex, *target_));
}

//----------------------------------------------------------------------------
int VectorReducer::copyIn(int numValues,
                          const int* indices, const double* values,
                          int vectorIndex)
{
  return(reducer_->addVectorValues(numValues, indices, values, false,
                                   isSolution_, vectorIndex, *target_));
}

//----------------------------------------------------------------------------
int VectorReducer::giveToUnderlyingVector(int numValues,
                                          const int* indices,
                                          const double* values,
                                          bool sumInto,
                                          int vectorIndex)
{
  int err = reducer_->addVectorValues(numValues, indices, values, sumInto,
                                      isSolution_, vectorIndex, *target_);
  return(err);
}

//----------------------------------------------------------------------------
int VectorReducer::sumInFieldData(int fieldID,
                                  int idType,
                                  int numIDs,
                                  const int* IDs,
                                  const double* data,
                                  int vectorIndex)
{
  fei::SharedPtr<fei::VectorSpace> vspace = target_->getVectorSpace();
  int fieldSize = vspace->getFieldSize(fieldID);
  int numIndices = numIDs*fieldSize;
  std::vector<int> indices(numIndices);
  int err = vspace->getGlobalIndices(numIDs, IDs, idType, fieldID, &indices[0]);
  if (err != 0) {
    throw std::runtime_error("fei::VectorReducer::sumInFieldData ERROR in vspace->getGlobalIndices.");
  }

  return(sumIn(numIndices, &indices[0], data, vectorIndex));
}

//----------------------------------------------------------------------------
int VectorReducer::copyInFieldData(int fieldID,
                                   int idType,
                                   int numIDs,
                                   const int* IDs,
                                   const double* data,
                                   int vectorIndex)
{
  fei::SharedPtr<fei::VectorSpace> vspace = target_->getVectorSpace();
  int fieldSize = vspace->getFieldSize(fieldID);
  int numIndices = numIDs*fieldSize;
  std::vector<int> indices(numIndices);
  int err = vspace->getGlobalIndices(numIDs, IDs, idType, fieldID, &indices[0]);
  if (err != 0) {
    throw std::runtime_error("fei::VectorReducer::copyInFieldData ERROR in vspace->getGlobalIndices.");
  }

  return(copyIn(numIndices, &indices[0], data, vectorIndex));
}

//----------------------------------------------------------------------------
int VectorReducer::copyInFieldDataLocalIDs(int fieldID,
                                   int idType,
                                   int numIDs,
                                   const int* localIDs,
                                   const double* data,
                                   int vectorIndex)
{
  fei::SharedPtr<fei::VectorSpace> vspace = target_->getVectorSpace();
  int fieldSize = vspace->getFieldSize(fieldID);
  int numIndices = numIDs*fieldSize;
  std::vector<int> indices(numIndices);
  int err = vspace->getGlobalIndicesLocalIDs(numIDs, localIDs, idType, fieldID, &indices[0]);
  if (err != 0) {
    throw std::runtime_error("fei::VectorReducer::copyInFieldData ERROR in vspace->getGlobalIndices.");
  }

  return(copyIn(numIndices, &indices[0], data, vectorIndex));
}

//----------------------------------------------------------------------------
int VectorReducer::copyOut_FE(int nodeNumber, int dofOffset,
                              double& value)
{
  return(-1);
}

//----------------------------------------------------------------------------
int VectorReducer::copyOutFieldData(int fieldID,
                                    int idType,
                                    int numIDs,
                                    const int* IDs,
                                    double* data,
                                    int vectorIndex)
{
  fei::SharedPtr<fei::VectorSpace> vspace = target_->getVectorSpace();
  int fieldSize = vspace->getFieldSize(fieldID);
  int numIndices = numIDs*fieldSize;
  std::vector<int> indices(numIndices);
  int err = vspace->getGlobalIndices(numIDs, IDs, idType, fieldID, &indices[0]);
  if (err != 0) {
    throw std::runtime_error("fei::VectorReducer::copyOutFieldData ERROR in vspace->getGlobalIndices.");
  }

  return(copyOut(numIndices, &indices[0], data, vectorIndex));
}

//----------------------------------------------------------------------------
int VectorReducer::writeToFile(const char* filename,
                               bool matrixMarketFormat)
{
  return( target_->writeToFile(filename, matrixMarketFormat) );
}

//----------------------------------------------------------------------------
int VectorReducer::writeToStream(FEI_OSTREAM& ostrm,
                                 bool matrixMarketFormat)
{
  return( target_->writeToStream(ostrm, matrixMarketFormat) );
}

//----------------------------------------------------------------------------
int VectorReducer::copyOut(int numValues,
                           const int* indices,
                           double* values,
                           int vectorIndex) const
{
  int err = reducer_->copyOutVectorValues(numValues, indices, values,
                                          isSolution_, vectorIndex,
                                          *target_);
  return(err);
}

//----------------------------------------------------------------------------
int VectorReducer::sumIntoFEVector(int blockID,
                                   int connOffset,
                                   int numNodes,
                                   const int* nodeNumbers,
                                   const int* numIndicesPerNode,
                                   const double* values)
{
  return(-1);
}

}//namespace fei

