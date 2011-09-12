/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "EpetraExt_MultiPointModelEvaluator.h"
#include "Epetra_Map.h"
#include "Teuchos_as.hpp"

EpetraExt::MultiPointModelEvaluator::MultiPointModelEvaluator(
    Teuchos::RefCountPtr<EpetraExt::ModelEvaluator> underlyingME_,
    const Teuchos::RefCountPtr<EpetraExt::MultiComm> &globalComm_,
    const std::vector<Epetra_Vector*> initGuessVec_,
    Teuchos::RefCountPtr<std::vector< Teuchos::RefCountPtr<Epetra_Vector> > >  q_vec_,
    Teuchos::RefCountPtr<std::vector< Teuchos::RefCountPtr<Epetra_Vector> > >  matching_vec_
    ) : 
    underlyingME(underlyingME_),
    globalComm(globalComm_),
    underlyingNg(0),
    timeStepsOnTimeDomain(globalComm_->NumTimeStepsOnDomain()),
    numTimeDomains(globalComm_->NumSubDomains()),
    timeDomain(globalComm_->SubDomainRank()),
    rowStencil(0),
    rowIndex(0),
    q_vec(q_vec_),
    matching_vec(matching_vec_)
{
  using Teuchos::as;
  if (globalComm->MyPID()==0) {
     cout  << "----------MultiPoint Partition Info------------"
           << "\n\tNumProcs              = " << globalComm->NumProc()
           << "\n\tSpatial Decomposition = " << globalComm->SubDomainComm().NumProc()
           << "\n\tNumber of Domains     = " << numTimeDomains
           << "\n\tSteps on Domain 0     = " << timeStepsOnTimeDomain
           << "\n\tTotal Number of Steps = " << globalComm->NumTimeSteps();
    cout   << "\n-----------------------------------------------" << endl;
    }

   // Construct global block matrix graph from split W and stencil,
   // which is just diagonal in this case

   split_W = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(underlyingME->create_W());

   rowStencil = new std::vector< std::vector<int> >(timeStepsOnTimeDomain);
   rowIndex = new std::vector<int>;
   for (int i=0; i < timeStepsOnTimeDomain; i++) {
     (*rowStencil)[i].push_back(0);
     (*rowIndex).push_back(i + globalComm->FirstTimeStepOnDomain());
   }

   block_W = Teuchos::rcp(new EpetraExt::BlockCrsMatrix(*split_W,
                               *rowStencil, *rowIndex, *globalComm));

   // Test for g vector
   EpetraExt::ModelEvaluator::OutArgs underlyingOutArgs = underlyingME->createOutArgs();

   underlyingNg = underlyingOutArgs.Ng();
   if (underlyingNg) {
     if (underlyingOutArgs.supports(OUT_ARG_DgDp,0,0).supports(DERIV_TRANS_MV_BY_ROW))
       orientation_DgDp = DERIV_TRANS_MV_BY_ROW;
     else
       orientation_DgDp = DERIV_MV_BY_COL;
   }

   // This code assumes 2 parameter vectors, 1 for opt, second for MultiPoint states
   TEST_FOR_EXCEPT(underlyingOutArgs.Np()!=2);

   // temporary quantities
   const Epetra_Map& split_map = split_W->RowMatrixRowMap();
   num_p0 =  underlyingME_->get_p_map(0)->NumMyElements();
   if (underlyingNg)  num_g0 = underlyingME_->get_g_map(0)->NumMyElements();
   else num_g0 = 0;
   num_dg0dp0 = num_g0 * num_p0;

   // Construct global solution vector, residual vector -- local storage
   block_x = new EpetraExt::BlockVector(split_map, block_W->RowMap());
   block_f = new EpetraExt::BlockVector(*block_x); 
   block_DfDp = new EpetraExt::BlockMultiVector(split_map, block_W->RowMap(), num_p0);
    if (underlyingNg)  
   block_DgDx = new EpetraExt::BlockMultiVector(split_map, block_W->RowMap(), num_g0);

   // Allocate local storage of epetra vectors
   split_x = Teuchos::rcp(new Epetra_Vector(split_map));
   split_f = Teuchos::rcp(new Epetra_Vector(split_map));
   split_DfDp = Teuchos::rcp(new Epetra_MultiVector(split_map, num_p0));
   if (underlyingNg)  
     split_DgDx = Teuchos::rcp(new Epetra_MultiVector(split_map, num_g0));
   if (underlyingNg) { 
     if(orientation_DgDp == DERIV_TRANS_MV_BY_ROW)
       split_DgDp = Teuchos::rcp(new Epetra_MultiVector(*(underlyingME_->get_p_map(0)), num_g0));
     else
       split_DgDp = Teuchos::rcp(new Epetra_MultiVector(*(underlyingME_->get_g_map(0)), num_p0));
   } 
   if (underlyingNg)  
     split_g = Teuchos::rcp(new Epetra_Vector(*(underlyingME_->get_g_map(0))));

   // Packaging required for getting multivectors back as Derivatives
   derivMV_DfDp = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DfDp);
   deriv_DfDp = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DfDp);
   if (underlyingNg)  {
     derivMV_DgDx = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DgDx, DERIV_TRANS_MV_BY_ROW);
     deriv_DgDx = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DgDx);
     derivMV_DgDp = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DgDp, orientation_DgDp);
     deriv_DgDp = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DgDp);
   }

   // For 4D, we will need the overlap vector and importer between them
   // Overlap not needed for MultiPoint -- no overlap between blocks
   /*   solutionOverlap = new EpetraExt::BlockVector(split_W->RowMatrixRowMap(),
                                                     block_W->ColMap());
        overlapImporter = new Epetra_Import(solutionOverlap->Map(), block_x->Map());
   */

   // Load initial guess into block solution vector
   solution_init = Teuchos::rcp(new EpetraExt::BlockVector(*block_x));
   for (int i=0; i < timeStepsOnTimeDomain; i++)
           solution_init->LoadBlockValues(*(initGuessVec_[i]), (*rowIndex)[i]);

 
   //Prepare logic for matching problem
   if (Teuchos::is_null(matching_vec))  matchingProblem = false;
   else matchingProblem = true;

   if (matchingProblem) {
     TEST_FOR_EXCEPT(as<int>(matching_vec->size())!=timeStepsOnTimeDomain);
     TEST_FOR_EXCEPT(!(*matching_vec)[0]->Map().SameAs(*(underlyingME_->get_g_map(0))));
     TEST_FOR_EXCEPT(num_g0 != 1); //This restriction may be lifted later
   }
}

EpetraExt::MultiPointModelEvaluator::~MultiPointModelEvaluator()
{
  delete block_x;
  delete block_f;
  delete block_DfDp;
  if (underlyingNg)  delete block_DgDx;
  delete rowStencil;
  delete rowIndex;

  delete derivMV_DfDp;
  delete deriv_DfDp;
  if (underlyingNg) {
    delete derivMV_DgDx;
    delete deriv_DgDx;
    delete derivMV_DgDp;
    delete deriv_DgDp;
  }
}

Teuchos::RefCountPtr<const Epetra_Map> EpetraExt::MultiPointModelEvaluator::get_x_map() const
{
  return Teuchos::rcp(&(block_W->OperatorDomainMap()), false);
}

Teuchos::RefCountPtr<const Epetra_Map> EpetraExt::MultiPointModelEvaluator::get_f_map() const
{
  return get_x_map();
}

Teuchos::RefCountPtr<const Epetra_Map> EpetraExt::MultiPointModelEvaluator::get_p_map(int l) const
{
  return underlyingME->get_p_map(l);
}

Teuchos::RefCountPtr<const Epetra_Map> EpetraExt::MultiPointModelEvaluator::get_g_map(int j) const
{
  return underlyingME->get_g_map(j);
}

Teuchos::RefCountPtr<const Epetra_Vector> EpetraExt::MultiPointModelEvaluator::get_x_init() const
{
  return solution_init;
}

Teuchos::RefCountPtr<const Epetra_Vector> EpetraExt::MultiPointModelEvaluator::get_p_init(int l) const
{
  return underlyingME->get_p_init(l);
}

Teuchos::RefCountPtr<Epetra_Operator> EpetraExt::MultiPointModelEvaluator::create_W() const
{
  return block_W;
}

EpetraExt::ModelEvaluator::InArgs EpetraExt::MultiPointModelEvaluator::createInArgs() const
{
  //return underlyingME->createInArgs();
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.set_Np(1);
  inArgs.setSupports(IN_ARG_x,true);
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs EpetraExt::MultiPointModelEvaluator::createOutArgs() const
{
  //return underlyingME->createOutArgs();
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, underlyingNg);
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_FULL
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DfDp,0,DERIV_MV_BY_COL);
  outArgs.set_DfDp_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_CONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );

  if (underlyingNg) {
    outArgs.setSupports(OUT_ARG_DgDx,0,DERIV_TRANS_MV_BY_ROW);
    outArgs.set_DgDx_properties(
      0,DerivativeProperties(
        DERIV_LINEARITY_NONCONST
        ,DERIV_RANK_DEFICIENT
        ,true // supportsAdjoint
        )
      );
    outArgs.setSupports(OUT_ARG_DgDp,0,0, orientation_DgDp);
    outArgs.set_DgDp_properties(
      0,0,DerivativeProperties(
        DERIV_LINEARITY_NONCONST
        ,DERIV_RANK_DEFICIENT
        ,true // supportsAdjoint
        )
      );
  }
  return outArgs;
}

void EpetraExt::MultiPointModelEvaluator::evalModel( const InArgs& inArgs,
                                            const OutArgs& outArgs ) const
{

  EpetraExt::ModelEvaluator::InArgs  underlyingInArgs  = underlyingME->createInArgs();
  EpetraExt::ModelEvaluator::OutArgs underlyingOutArgs = underlyingME->createOutArgs();

  //temp code for multipoint param q vec
/*
  Teuchos::RefCountPtr<Epetra_Vector> q =
    Teuchos::rcp(new Epetra_Vector(*(underlyingME->get_p_map(1))));
*/

  // Parse InArgs
  Teuchos::RefCountPtr<const Epetra_Vector> p_in = inArgs.get_p(0);
  if (p_in.get()) underlyingInArgs.set_p(0, p_in);

  Teuchos::RefCountPtr<const Epetra_Vector> x_in = inArgs.get_x();
  block_x->Epetra_Vector::operator=(*x_in); //copy into block vector

  // Parse OutArgs
  Teuchos::RefCountPtr<Epetra_Vector> f_out = outArgs.get_f();

  Teuchos::RefCountPtr<Epetra_Operator> W_out = outArgs.get_W();
  Teuchos::RefCountPtr<EpetraExt::BlockCrsMatrix> W_block =
     Teuchos::rcp_dynamic_cast<EpetraExt::BlockCrsMatrix>(W_out);

  Teuchos::RefCountPtr<Epetra_Vector> g_out;
  if (underlyingNg) g_out = outArgs.get_g(0); 
  if (g_out.get()) g_out->PutScalar(0.0);

  EpetraExt::ModelEvaluator::Derivative DfDp_out = outArgs.get_DfDp(0);

  EpetraExt::ModelEvaluator::Derivative DgDx_out;
  EpetraExt::ModelEvaluator::Derivative DgDp_out;
  if (underlyingNg) {
    DgDx_out = outArgs.get_DgDx(0);
    DgDp_out = outArgs.get_DgDp(0,0);
    if (!DgDx_out.isEmpty()) DgDx_out.getMultiVector()->PutScalar(0.0);
    if (!DgDp_out.isEmpty()) DgDp_out.getMultiVector()->PutScalar(0.0);
  }

  // For mathcingProblems, g is needed to calc DgDx DgDp, so ask for
  //  g even if it isn't requested.
  bool need_g = g_out.get();
  if (matchingProblem)
    if ( !DgDx_out.isEmpty() || !DgDp_out.isEmpty() ) need_g = true;


  // Begin loop over Points (steps) owned on this proc
  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    // Set MultiPoint parameter vector
    underlyingInArgs.set_p(1, (*q_vec)[i]);

    // Set InArgs
    block_x->ExtractBlockValues(*split_x, (*rowIndex)[i]);
    underlyingInArgs.set_x(split_x);

    // Set OutArgs
    if (f_out.get()) underlyingOutArgs.set_f(split_f);

    if (need_g) underlyingOutArgs.set_g(0, split_g);

    if (W_out.get()) underlyingOutArgs.set_W(split_W);

    if (!DfDp_out.isEmpty()) underlyingOutArgs.set_DfDp(0, *deriv_DfDp);

    if (!DgDx_out.isEmpty()) underlyingOutArgs.set_DgDx(0, *deriv_DgDx);
  
    if (!DgDp_out.isEmpty()) underlyingOutArgs.set_DgDp(0, 0, *deriv_DgDp);

    //********Eval Model ********/
    underlyingME->evalModel(underlyingInArgs, underlyingOutArgs);
    //********Eval Model ********/

    // If matchingProblem, modify all g-related quantitites G = 0.5*(g-g*)^2 / g*^2
    if (matchingProblem) {
      if (need_g) {
        double diff = (*split_g)[0] -  (*(*matching_vec)[i])[0];
        double nrmlz = fabs((*(*matching_vec)[i])[0]) + 1.0e-6;
        (*split_g)[0] = 0.5 * diff * diff/(nrmlz*nrmlz);
        if (!DgDx_out.isEmpty()) split_DgDx->Scale(diff/(nrmlz*nrmlz));
        if (!DgDp_out.isEmpty()) split_DgDp->Scale(diff/(nrmlz*nrmlz));
      }
    }

    // Repackage block components into global block matrx/vector/multivector
    if (f_out.get()) block_f->LoadBlockValues(*split_f, (*rowIndex)[i]);
    if (W_out.get()) W_block->LoadBlock(*split_W, i, 0);
        // note: split_DfDp points inside deriv_DfDp
    if (!DfDp_out.isEmpty()) block_DfDp->LoadBlockValues(*split_DfDp, (*rowIndex)[i]);
    if (!DgDx_out.isEmpty()) block_DgDx->LoadBlockValues(*split_DgDx, (*rowIndex)[i]);

    // Assemble multiple steps on this domain into g and dgdp(0) vectors
    if (g_out.get()) g_out->Update(1.0, *split_g, 1.0);

    if (!DgDp_out.isEmpty())
      DgDp_out.getMultiVector()->Update(1.0, *split_DgDp, 1.0);

  } // End loop over multiPoint steps on this domain/cluster

  //Copy block vectors into *_out vectors of same size
  if (f_out.get()) f_out->operator=(*block_f);
  if (!DfDp_out.isEmpty()) 
    DfDp_out.getMultiVector()->operator=(*block_DfDp);
  if (!DgDx_out.isEmpty()) 
    DgDx_out.getMultiVector()->operator=(*block_DgDx);

  //Sum together obj fn contributions from differnt Domains (clusters).
  if (numTimeDomains > 1) {
    double factorToZeroOutCopies = 0.0;
    if (globalComm->SubDomainComm().MyPID()==0) factorToZeroOutCopies = 1.0;
    if (g_out.get()) {
      (*g_out).Scale(factorToZeroOutCopies);
      double* vPtr = &((*g_out)[0]);
      Epetra_Vector tmp = *(g_out.get());
      globalComm->SumAll( &(tmp[0]), vPtr, num_g0);
    }
    if (!DgDp_out.isEmpty()) {
      DgDp_out.getMultiVector()->Scale(factorToZeroOutCopies);
      double* mvPtr = (*DgDp_out.getMultiVector())[0];
      Epetra_MultiVector tmp = *(DgDp_out.getMultiVector());
      globalComm->SumAll(tmp[0], mvPtr, num_dg0dp0);
    }
  }
}
