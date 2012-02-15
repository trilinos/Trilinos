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


#include <fei_macros.hpp>

#include <fei_MatrixGraph.hpp>
#include <fei_LinearSystem.hpp>
#include <snl_fei_LinearSystem_General.hpp>
#include <snl_fei_Utils.hpp>

//----------------------------------------------------------------------------
fei::LinearSystem::LinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
 : matrix_(),
   soln_(),
   rhs_(),
   matrixGraph_(matrixGraph),
   dbcManager_(NULL)
{
}

//----------------------------------------------------------------------------
fei::LinearSystem::~LinearSystem()
{
  delete dbcManager_;

  for(unsigned i=0; i<attributeNames_.size(); ++i) {
    delete [] attributeNames_[i];
  }
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::LinearSystem>
fei::LinearSystem::Factory::createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
{
  fei::SharedPtr<fei::LinearSystem>
    linsys(new snl_fei::LinearSystem_General(matrixGraph));

  return(linsys);
}

//----------------------------------------------------------------------------
void fei::LinearSystem::setMatrix(fei::SharedPtr<fei::Matrix>& matrix)
{
  matrix_ = matrix;
}

//----------------------------------------------------------------------------
int fei::LinearSystem::putAttribute(const char* name,
                                    void* attribute)
{
  snl_fei::storeNamedAttribute(name, attribute,
                               attributeNames_, attributes_);
  return(0);
}

//----------------------------------------------------------------------------
int fei::LinearSystem::getAttribute(const char* name,
                                    void*& attribute)
{
  attribute = snl_fei::retrieveNamedAttribute(name, attributeNames_, attributes_);
  return(attribute==NULL ? -1 : 0);
}

//----------------------------------------------------------------------------
int fei::LinearSystem::loadEssentialBCs(int numIDs,
                                 const int* IDs,
                                 int idType,
                                 int fieldID,
                                 int offsetIntoField,
                                 const double* prescribedValues)
{
  if (dbcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager(matrixGraph_->getRowSpace());
  }

  try {
    dbcManager_->addBCRecords(numIDs, idType, fieldID, offsetIntoField,
                              IDs, prescribedValues);
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what()<<FEI_ENDL;
    return(-1);
  }

  return(0);
}

//----------------------------------------------------------------------------
int fei::LinearSystem::loadEssentialBCs(int numIDs,
                                 const int* IDs,
                                 int idType,
                                 int fieldID,
                                 const int* offsetsIntoField,
                                 const double* prescribedValues)
{
  if (dbcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager(matrixGraph_->getRowSpace());
  }

  try {
    dbcManager_->addBCRecords(numIDs, idType, fieldID, IDs, offsetsIntoField,
                              prescribedValues);
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what()<<FEI_ENDL;
    return(-1);
  }

  return(0);
}

