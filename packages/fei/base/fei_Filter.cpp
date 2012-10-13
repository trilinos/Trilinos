/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_sstream.hpp"

#include "fei_CommUtils.hpp"
#include "fei_TemplateUtils.hpp"

#include "fei_defs.h"
#include "fei_NodeDescriptor.hpp"
#include "fei_NodeDatabase.hpp"
#include "fei_BlockDescriptor.hpp"
#include "SNL_FEI_Structure.hpp"
#include "snl_fei_Utils.hpp"
#include "fei_Filter.hpp"

#include <cmath>
#include <algorithm>

#undef fei_file
#define fei_file "fei_Filter.cpp"

#include "fei_ErrMacros.hpp"

//------------------------------------------------------------------------------
Filter::Filter(SNL_FEI_Structure* probStruct)
  : problemStructure_(probStruct),
    logInput_(false),
    logInputStream_(NULL),
    outputLevel_(0)
{
}

//------------------------------------------------------------------------------
Filter::~Filter()
{}

//------------------------------------------------------------------------------
void Filter::setLogStream(FEI_OSTREAM* logstrm)
{
  logInputStream_ = logstrm;
}

//------------------------------------------------------------------------------
FEI_OSTREAM* Filter::logStream()
{
  return( logInputStream_ );
}

//------------------------------------------------------------------------------
void Filter::copyStiffness(const double* const* elemStiff,
			   int numRows, int elemFormat,
			   double** copy)
{
  //
  //Unpack the element stiffness array in elemStiff into a full dense
  //'copy'.
  //
  int i, j;
  const double* elStiff_i = NULL;

  switch (elemFormat) {

  case FEI_DENSE_ROW:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      for (j = 0; j < numRows; j++) {
	copy[i][j] = elStiff_i[j];
      }
    }
    break;

  case FEI_UPPER_SYMM_ROW:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      int jcol=0;
      for (j = i; j < numRows; j++) {
	copy[i][j] = elStiff_i[jcol++];
	copy[j][i] = copy[i][j];
      }
    }
    break;

  case FEI_LOWER_SYMM_ROW:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      for (j = 0; j <=i; j++) {
	copy[i][j] = elStiff_i[j];
	copy[j][i] = copy[i][j];
      }
    }
    break;

  case FEI_DENSE_COL:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      for (j = 0; j < numRows; j++) {
	copy[j][i] = elStiff_i[j];
      }
    }
    break;

  case FEI_UPPER_SYMM_COL:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      for (j = 0; j <= i; j++) {
	copy[i][j] = elStiff_i[j];
	copy[j][i] = copy[i][j];
      }
    }
    break;

  case FEI_LOWER_SYMM_COL:
    for (i = 0; i < numRows; i++) {
      elStiff_i = elemStiff[i];
      int jcol=0;
      for (j = i; j < numRows; j++) {
	copy[i][j] = elStiff_i[jcol++];
	copy[j][i] = copy[i][j];
      }
    }
    break;

  default:
    throw std::runtime_error("copyStiffness ERROR, unrecognized elem-format");
  }
}

//------------------------------------------------------------------------------
const NodeDescriptor* Filter::findNode(GlobalID nodeID) const {
//
//This function returns a NodeDescriptor ptr, may return NULL.
//
  const NodeDescriptor* node = NULL;
  problemStructure_->getNodeDatabase().getNodeWithID(nodeID, node);
  return node;
}

//------------------------------------------------------------------------------
const NodeDescriptor& Filter::findNodeDescriptor(GlobalID nodeID) const {
//
//This function returns a NodeDescriptor reference if nodeID is an active node.
//
  const NodeDescriptor* node = NULL;
  int err = problemStructure_->getNodeDatabase().getNodeWithID(nodeID, node);

  if (err != 0) {
    fei::console_out() << "ERROR, Filter::findNodeDescriptor unable to find node "
	 << static_cast<int>(nodeID) << FEI_ENDL;
    std::abort();
  }

  return( *node );
}
//------------------------------------------------------------------------------
int Filter::calculateResidualNorms(int whichNorm, int numFields,
				   int* fieldIDs, double* norms,
				   std::vector<double>& residValues)
{
  std::vector<double> normsArray(numFields, 0.0);

  std::fill(fieldIDs, fieldIDs+numFields, -999);

  std::vector<double> tmpNorms(numFields);
  double* tmpNormsPtr = &tmpNorms[0];

  double* residPtr = &(residValues[0]);

  const std::vector<int>& pfieldIDs = problemStructure_->getFieldIDs();
  int numDBFields = pfieldIDs.size();
  std::vector<int>::const_iterator
    f_iter = pfieldIDs.begin(),
    f_end = pfieldIDs.end();

  int DBFieldSize = 0;

  int offset = 0;
  for(; f_iter != f_end; ++f_iter) {

    if (offset == 0) DBFieldSize = problemStructure_->getFieldSize(*f_iter);

    if (*f_iter > -1) {
      if (offset < numFields) {
        fieldIDs[offset] = *f_iter;
        tmpNormsPtr[offset++] = 0.0;
      }
    }
  }

  int reducedStartRow = problemStructure_->getFirstReducedEqn();
  int reducedEndRow   = problemStructure_->getLastReducedEqn();

  NodeDatabase& nodeDB = problemStructure_->getNodeDatabase();
  int numNodes = nodeDB.getNumNodeDescriptors();

  bool haveSlaves = problemStructure_->numSlaveEquations() > 0;

  for(int i=0; i<numNodes; i++) {
    const NodeDescriptor* node = NULL;
    nodeDB.getNodeAtIndex(i, node);

    if (node==NULL || node->getOwnerProc() != localRank_) continue;

    const int* fieldIDList = node->getFieldIDList();
    const int* fieldEqnNums = node->getFieldEqnNumbers();
    int numNodeFields = node->getNumFields();

    for(int j=0; j<numNodeFields; j++) {
      int fIndex = 0;
      int fSize = DBFieldSize;

      if (numDBFields > 1) {
        fIndex = fei::binarySearch(fieldIDList[j], fieldIDs, numFields);
        if (fIndex < 0) return(-1);
        fSize = problemStructure_->getFieldSize(fieldIDList[j]);
        if (fSize < 0) return(-1);
      }

      for(int k=0; k<fSize; k++) {
        int eqn = fieldEqnNums[j]+k;

        if (haveSlaves) {
          if (problemStructure_->isSlaveEqn(eqn)) continue;
          int reducedEqn;
          problemStructure_->translateToReducedEqn(eqn, reducedEqn);

          if (reducedEqn < reducedStartRow || reducedEqn > reducedEndRow) {
            continue;
          }
          eqn = reducedEqn;
        }

        double rval = residPtr[eqn - reducedStartRow];

        switch(whichNorm) {
          case 0:
            if (tmpNormsPtr[fIndex] < std::abs(rval)) tmpNormsPtr[fIndex] = rval;
            break;
          case 1:
            tmpNormsPtr[fIndex] += std::abs(rval);
            break;
          case 2:
            tmpNormsPtr[fIndex] += rval*rval;
            break;
          default:
            FEI_COUT << "Filter::residualNorm: ERROR, whichNorm="<<whichNorm
              << " not recognized." << FEI_ENDL;
            return(-1);
        }
      }
    }
  }

  //so at this point we have the local norms. We now need to perform the
  //global max or sum, depending on which norm it is.

  MPI_Comm comm = problemStructure_->getCommunicator();

  if (whichNorm != 0) {
    CHK_ERR( fei::GlobalSum(comm, tmpNorms, normsArray) );
  }
  else {
    CHK_ERR( fei::GlobalMax(comm, tmpNorms, normsArray) );
  }

  for(int i=0; i<numFields; ++i) {
    norms[i] = normsArray[i];
  }

  if (whichNorm == 2) {
    for(int i=0; i<numFields; ++i) norms[i] = std::sqrt(norms[i]);
  }

  return(0);
}

//------------------------------------------------------------------------------
int Filter::parameters(int numParams, const char *const* paramStrings)
{
  if (numParams == 0 || paramStrings == NULL) return(0);

  const char* param = snl_fei::getParam("outputLevel",numParams,paramStrings);

  if ( param != NULL){
    std::string str(&(param[11]));
    FEI_ISTRINGSTREAM isstr(str);
    isstr >> outputLevel_;
  }

  return(0);
}

