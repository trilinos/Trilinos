/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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

int
DirichletBCManager::getEqnNumber(int IDType, int ID, int fieldID, int offsetIntoField)
{
  int eqn = -1;
  try {
    if (vecSpace_.get() != NULL) {
      vecSpace_->getGlobalIndex(IDType, ID, fieldID, eqn);
    }
    else {
      if (structure_ == NULL) {
        throw std::runtime_error("fei::DirichletBCManager has NULL SNL_FEI_Structure.");
      }
      NodeDatabase& nodeDB = structure_->getNodeDatabase();
      const NodeDescriptor* node = NULL;
      nodeDB.getNodeWithID(ID, node);
      if (node == NULL) {
        throw std::runtime_error("fei::DirichletBCManager::getEqnNumber failed to get node.");
      }
      node->getFieldEqnNumber(fieldID, eqn);
    }
  }
  catch(std::runtime_error& exc) {
    FEI_OSTRINGSTREAM osstr;
    osstr << "fei::DirichletBCManager::finalizeBCEqns caught exception: "
       << exc.what() << " BC IDType="<<IDType<<", ID="<<ID
       << ", fieldID="<<fieldID;
    fei::console_out() << osstr.str() << FEI_ENDL;
    ERReturn(-1);
  }

  return eqn + offsetIntoField;
}

void
DirichletBCManager::addBCRecords(int numBCs,
                                 int IDType,
                                 int fieldID,
                                 int offsetIntoField,
                                 const int* IDs,
                                 const double* prescribedValues)
{
  for(int i=0; i<numBCs; ++i) {
    int eqn = getEqnNumber(IDType, IDs[i], fieldID, offsetIntoField);

    bc_map::iterator iter = bcs_.lower_bound(eqn);

    if (iter == bcs_.end() || iter->first != eqn) {
      bcs_.insert(iter, std::make_pair(eqn, prescribedValues[i]));
      continue;
    }
    else iter->second = prescribedValues[i];
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
    int eqn = getEqnNumber(IDType, IDs[i], fieldID, offsetsIntoField[i]);

    bc_map::iterator iter = bcs_.lower_bound(eqn);

    if (iter == bcs_.end() || iter->first != eqn) {
      bcs_.insert(iter, std::make_pair(eqn, prescribedValues[i]));
      continue;
    }
    else iter->second = prescribedValues[i];
  }
}

int
DirichletBCManager::finalizeBCEqns(fei::Matrix& matrix,
                                   bool throw_if_bc_slave_conflict)
{
  fei::SharedPtr<fei::Reducer> reducer = matrix.getMatrixGraph()->getReducer();
  bool haveSlaves = reducer.get()!=NULL;

  //copy the boundary-condition prescribed values into the matrix, in
  //an equation-number obtained by using the matrix' VectorSpace to map
  //from the BC's idtype,id,fieldID,component to an equation-number. The
  //bc values will go on the diagonal of the matrix, i.e., column-index
  //will be the same equation-number.

  bc_map::iterator iter = bcs_.begin(), iter_end = bcs_.end();

  for(; iter!=iter_end; ++iter) {

    int eqn = iter->first;

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

    double* ptr = &iter->second;

    CHK_ERR( matrix.copyIn(1, &eqn, 1, &eqn, &ptr) );
  }

  bcs_.clear();
  return(0);
}

int
DirichletBCManager::finalizeBCEqns(EqnBuffer& bcEqns)
{
  //copy the boundary-condition prescribed values into bcEqns.

  bc_map::iterator iter = bcs_.begin(), iter_end = bcs_.end();

  for(; iter!=iter_end; ++iter) {
    int eqn = iter->first;
    double coef = iter->second;

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

