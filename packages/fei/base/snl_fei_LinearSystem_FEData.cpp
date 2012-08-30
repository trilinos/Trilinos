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

#include <math.h>

#include <fei_macros.hpp>

#include <fei_FiniteElementData.hpp>

#include <fei_CommUtils.hpp>
#include <snl_fei_LinearSystem_FEData.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix_Impl.hpp>
#include <snl_fei_Constraint.hpp>
#include <snl_fei_Utils.hpp>
#include <fei_impl_utils.hpp>

#include <fei_DirichletBCRecord.hpp>
#include <fei_DirichletBCManager.hpp>
#include <fei_EqnBuffer.hpp>
#include <fei_FEDataFilter.hpp>

#undef fei_file
#define fei_file "snl_fei_LinearSystem_FEData.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
snl_fei::LinearSystem_FEData::LinearSystem_FEData(fei::SharedPtr<FiniteElementData>& feData,
						  fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
  : fei::LinearSystem(matrixGraph),
    comm_(matrixGraph->getRowSpace()->getCommunicator()),
    localProc_(0),
    numProcs_(1),
    feData_(feData),
    lookup_(NULL)
{
  localProc_ = fei::localProc(comm_);
  numProcs_  = fei::numProcs(comm_);
}

//----------------------------------------------------------------------------
snl_fei::LinearSystem_FEData::~LinearSystem_FEData()
{
}

//----------------------------------------------------------------------------
bool snl_fei::LinearSystem_FEData::eqnIsEssentialBC(int globalEqnIndex) const
{
  fei::console_out() << "LinearSystem_FEData::eqnIsEssentialBC NOT IMPLEMENTED!!"<<FEI_ENDL;
  return(false);
}

//----------------------------------------------------------------------------
void snl_fei::LinearSystem_FEData::getEssentialBCs(std::vector<int>& bcEqns,
                                              std::vector<double>& bcVals) const
{
  fei::console_out() << "LinearSystem_FEData::getEssentialBC_Eqns NOT IMPLEMENTED!!"<<FEI_ENDL;
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
  if (dbcManager_ == NULL) {
    dbcManager_ = new fei::DirichletBCManager(matrixGraph_->getRowSpace());
  }

  int err = 0;
  if (matrix_.get() != NULL && globalAssemble) {
    err = matrix_->gatherFromOverlap();
    if (err != 0) {
      fei::console_out() << "snl_fei::LinearSystem_FEData::loadComplete, ERROR in matrix."
	   << "gatherFromOverlap(), data may be incorrect."<<FEI_ENDL;
    }
  }

  fei::Barrier(comm_);

  if (rhs_.get() != NULL && globalAssemble) {
    err = rhs_->gatherFromOverlap();
    if (err != 0) {
      fei::console_out() << "snl_fei::LinearSystem_FEData::loadComplete, ERROR rhs."
        << "gatherFromOverlap(), data may be incorrect."<<FEI_ENDL;
    }
  }

  fei::Barrier(comm_);

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
  std::vector<int> essEqns;
  std::vector<double> essGamma;

  fei::SharedPtr<fei::FillableMat> localBCeqns(new fei::FillableMat);
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph_->getRowSpace();
  int localsize = vecSpace->getNumIndices_Owned();
  bool zeroSharedRows = false;
  fei::Matrix_Impl<fei::FillableMat> bcEqns(localBCeqns, matrixGraph_, localsize, zeroSharedRows);

  CHK_ERR( dbcManager_->finalizeBCEqns(bcEqns) );

  std::map<int,fei::FillableMat*>& remotes = bcEqns.getRemotelyOwnedMatrices();
  std::map<int,fei::FillableMat*>::iterator
    it = remotes.begin(),
    it_end = remotes.end();
  for(; it!=it_end; ++it) {
    fei::impl_utils::separate_BC_eqns( *(it->second), essEqns, essGamma);
  }

  CHK_ERR( bcEqns.gatherFromOverlap() );

  fei::impl_utils::separate_BC_eqns( *(bcEqns.getMatrix()), essEqns, essGamma);

  int len = essEqns.size();
  essEqns.resize(len*3);
  int* nodeNumbers = &essEqns[0]+len;
  int* dof_ids = nodeNumbers+len;

  if (len > 0) {
    int* essEqnsPtr = &essEqns[0];
    double* gammaPtr = &essGamma[0];

    for(int i=0; i<len; ++i) {
      int eqn = essEqnsPtr[i];
      nodeNumbers[i] = lookup_->getAssociatedNodeNumber(eqn);
      int fieldID = lookup_->getAssociatedFieldID(eqn);
      int base_eqn = lookup_->getEqnNumber(nodeNumbers[i], fieldID);
      dof_ids[i]=vecSpace->getFieldDofMap().get_dof_id(fieldID, eqn-base_eqn);
    }

    if (applyBCs) {
      CHK_ERR( feData_->setDirichletBCs(len, nodeNumbers, dof_ids, gammaPtr) );
    }
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
