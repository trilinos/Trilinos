/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iostream.hpp>
#include <fei_sstream.hpp>
#include <fei_DirichletBCManager.hpp>
#include <fei_DirichletBCRecord.hpp>
#include <fei_NodeDatabase.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Matrix.hpp>

#include <algorithm>
#include <vector>

typedef std::vector<fei::DirichletBCRecord> DBCvec;

#undef fei_file
#define fei_file "fei_DirichletBCManager.cpp"
#include <fei_ErrMacros.hpp>

namespace fei {
void
DirichletBCManager::addBCRecords(int numBCs,
                                 int IDType,
                                 int fieldID,
                                 int offsetIntoField,
                                 const int* IDs,
                                 const double* prescribedValues)
{
  for(int i=0; i<numBCs; ++i) {
    DirichletBCRecord dbc;
    dbc.IDType = IDType;
    dbc.ID = IDs[i];
    dbc.fieldID = fieldID;
    dbc.whichComponent = offsetIntoField;
    dbc.prescribedValue = prescribedValues[i];

    DBCvec::iterator iter = std::lower_bound(bcs_.begin(), bcs_.end(),
                                             dbc,less_DirichletBCRecord());

    if (iter == bcs_.end() || *iter != dbc) {
      bcs_.insert(iter, dbc);
      continue;
    }
    else iter->prescribedValue = dbc.prescribedValue;
  }
}

void
DirichletBCManager::addBCRecords(int numBCs,
                                 int IDType,
                                 int fieldID,
                                 const int* IDs,
                                 const int* offsetsIntoField,
                                 const double* prescribedValues)
{
  for(int i=0; i<numBCs; ++i) {
    DirichletBCRecord dbc;
    dbc.IDType = IDType;
    dbc.ID = IDs[i];
    dbc.fieldID = fieldID;
    dbc.whichComponent = offsetsIntoField[i];
    dbc.prescribedValue = prescribedValues[i];

    DBCvec::iterator iter = std::lower_bound(bcs_.begin(), bcs_.end(),
                                             dbc,less_DirichletBCRecord());
    if (iter == bcs_.end() || *iter != dbc) {
      bcs_.insert(iter, dbc);
      continue;
    }
    else iter->prescribedValue = dbc.prescribedValue;
  }
}

int
DirichletBCManager::finalizeBCEqns(fei::Matrix& matrix,
                                   bool throw_if_bc_slave_conflict)
{
  fei::VectorSpace& vecSpace = *(matrix.getMatrixGraph()->getRowSpace());
  fei::SharedPtr<fei::Reducer> reducer = matrix.getMatrixGraph()->getReducer();
  bool haveSlaves = reducer.get()!=NULL;

  //copy the boundary-condition prescribed values into the matrix, in
  //an equation-number obtained by using the matrix' VectorSpace to map
  //from the BC's idtype,id,fieldID,component to an equation-number. The
  //bc values will go on the diagonal of the matrix, i.e., column-index
  //will be the same equation-number.

  for(size_t i=0; i<bcs_.size(); ++i) {
    DirichletBCRecord& dbc = bcs_[i];
    int eqn = -1;
    try {
      CHK_ERR(vecSpace.getGlobalIndex(dbc.IDType, dbc.ID, dbc.fieldID, eqn));
    }
    catch(std::runtime_error& exc) {
      FEI_OSTRINGSTREAM osstr;
      osstr << "fei::DirichletBCManager::finalizeBCEqns caught exception: "
        << exc.what() << " BC IDType="<<dbc.IDType<<", ID="<<dbc.ID
       << ", fieldID="<<dbc.fieldID;
      FEI_CERR << osstr.str() << FEI_ENDL;
      ERReturn(-1);
    }
    eqn += dbc.whichComponent;

    if (haveSlaves) {
      if (reducer->isSlaveEqn(eqn)) {
        if (throw_if_bc_slave_conflict) {
          FEI_OSTRINGSTREAM osstr;
          osstr << "fei BCManager::finalizeBCeqns ERROR, eqn="<<eqn
            << " is both a BC eqn and slave-constraint eqn.";
          throw std::runtime_error(osstr.str());
        }
        continue;
      }
    }

    double* ptr = &dbc.prescribedValue;

    CHK_ERR( matrix.copyIn(1, &eqn, 1, &eqn, &ptr) );
  }

  bcs_.clear();
  return(0);
}

int
DirichletBCManager::finalizeBCEqns(NodeDatabase& nodeDB,
                                   EqnBuffer& bcEqns)
{
  //copy the boundary-condition prescribed values into bcEqns.
  double coef = 0.0;

  for(unsigned i=0; i<bcs_.size(); ++i) {
    DirichletBCRecord& dbc = bcs_[i];
    NodeDescriptor* node = NULL;
    nodeDB.getNodeWithID(dbc.ID, node);
    int fieldID = dbc.fieldID;

    int eqn = 0;
    if (!node->getFieldEqnNumber(fieldID, eqn)) {
      ERReturn(-1);
    }

    eqn += dbc.whichComponent;
    coef = dbc.prescribedValue;

    CHK_ERR( bcEqns.addEqn(eqn, &coef, &eqn, 1, false) );
  }

  bcs_.clear();
  return(0);
}

size_t
DirichletBCManager::getNumBCRecords() const
{
  return bcs_.size();
}

void
DirichletBCManager::clearAllBCs()
{
  bcs_.clear();
}

}//namespace fei

