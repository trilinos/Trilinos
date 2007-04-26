#include "EpetraExt_MultiPointModelEvaluator.h"

EpetraExt::MultiPointModelEvaluator::MultiPointModelEvaluator(
    Teuchos::RefCountPtr<EpetraExt::ModelEvaluator> underlyingME_,
    const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> &globalComm_,
    const std::vector<Epetra_Vector*> initGuessVec_,
    Teuchos::RefCountPtr<std::vector< Teuchos::RefCountPtr<Epetra_Vector> > >  q_vec_
    ) : 
    underlyingME(underlyingME_),
    globalComm(globalComm_),
    timeStepsOnTimeDomain(globalComm_->NumTimeStepsOnDomain()),
    numTimeDomains(globalComm_->NumSubDomains()),
    timeDomain(globalComm_->SubDomainRank()),
    rowStencil(0),
    rowIndex(0),
    q_vec(q_vec_)
{
  if (globalComm->MyPID()==0) {
     // TODO: pass in some output stream
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

   // temporary quantities
   const Epetra_Map& split_map = split_W->RowMatrixRowMap();
   int  num_p0 =  underlyingME_->get_p_map(0)->NumMyElements();
   int  num_g0 =  underlyingME_->get_g_map(0)->NumMyElements();

   // Construct global solution vector, residual vector -- local storage
   block_x = new EpetraExt::BlockVector(split_map, block_W->RowMap());
   block_f = new EpetraExt::BlockVector(*block_x); 
   block_DfDp = new EpetraExt::BlockMultiVector(split_map, block_W->RowMap(), num_p0);
   block_DgDx = new EpetraExt::BlockMultiVector(split_map, block_W->RowMap(), num_g0);

   // Allocate local storage of epetra vectors
   split_x = Teuchos::rcp(new Epetra_Vector(split_map));
   split_f = Teuchos::rcp(new Epetra_Vector(split_map));
   split_DfDp = Teuchos::rcp(new Epetra_MultiVector(split_map, num_p0));
   split_DgDx = Teuchos::rcp(new Epetra_MultiVector(split_map, num_g0));
   split_DgDp = Teuchos::rcp(new Epetra_MultiVector(*(underlyingME_->get_p_map(0)), num_g0));
   split_g = Teuchos::rcp(new Epetra_Vector(*(underlyingME_->get_g_map(0))));

   // Packaging required for getting multivectors back as Derivatives
   derivMV_DfDp = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DfDp);
   deriv_DfDp = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DfDp);
   derivMV_DgDx = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DgDx, DERIV_TRANS_MV_BY_ROW);
   deriv_DgDx = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DgDx);
   derivMV_DgDp = new EpetraExt::ModelEvaluator::DerivativeMultiVector(split_DgDp, DERIV_TRANS_MV_BY_ROW);
   deriv_DgDp = new EpetraExt::ModelEvaluator::Derivative(*derivMV_DgDp);

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

}

EpetraExt::MultiPointModelEvaluator::~MultiPointModelEvaluator()
{
  delete block_x;
  delete block_f;
  delete rowStencil;
  delete rowIndex;

  delete derivMV_DfDp;
  delete deriv_DfDp;
  delete derivMV_DgDx;
  delete deriv_DgDx;
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
  outArgs.set_Np_Ng(1,1);
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
  outArgs.setSupports(OUT_ARG_DgDx,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDx_properties(
    0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
  outArgs.setSupports(OUT_ARG_DgDp,0,0,DERIV_TRANS_MV_BY_ROW);
  outArgs.set_DgDp_properties(
    0,0,DerivativeProperties(
      DERIV_LINEARITY_NONCONST
      ,DERIV_RANK_DEFICIENT
      ,true // supportsAdjoint
      )
    );
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

  Teuchos::RefCountPtr<Epetra_Vector> g_out = outArgs.get_g(0);
  if (g_out.get()) g_out->PutScalar(0.0);

  EpetraExt::ModelEvaluator::Derivative DfDp_out = outArgs.get_DfDp(0);
  EpetraExt::ModelEvaluator::Derivative DgDx_out = outArgs.get_DgDx(0);
  EpetraExt::ModelEvaluator::Derivative DgDp_out = outArgs.get_DgDp(0,0);
  if (!DgDp_out.isEmpty()) DgDp_out.getMultiVector()->PutScalar(0.0);

  // Begin loop over steps owned on this proc
  for (int i=0; i < timeStepsOnTimeDomain; i++) {

    // Set MultiPoint parameter vector
    underlyingInArgs.set_p(1, (*q_vec)[i]);

    // Set InArgs
    block_x->ExtractBlockValues(*split_x, (*rowIndex)[i]);
    underlyingInArgs.set_x(split_x);

    // Set OutArgs
    if (f_out.get()) underlyingOutArgs.set_f(split_f);

    if (g_out.get()) underlyingOutArgs.set_g(0, split_g);

    if (W_out.get()) underlyingOutArgs.set_W(split_W);

    if (!DfDp_out.isEmpty()) underlyingOutArgs.set_DfDp(0, *deriv_DfDp);

    if (!DgDx_out.isEmpty()) underlyingOutArgs.set_DgDx(0, *deriv_DgDx);
  
    if (!DgDp_out.isEmpty()) underlyingOutArgs.set_DgDp(0, 0, *deriv_DgDp);

    //********Eval Model ********/
    underlyingME->evalModel(underlyingInArgs, underlyingOutArgs);
    //********Eval Model ********/

    // Repackage block components into global block matrx/vector/multivector
    if (f_out.get()) block_f->LoadBlockValues(*split_f, (*rowIndex)[i]);
    if (W_out.get()) W_block->LoadBlock(*split_W, i, 0);
        // note: split_DfDp points inside deriv_DfDp
    if (!DfDp_out.isEmpty()) block_DfDp->LoadBlockValues(*split_DfDp, (*rowIndex)[i]);
    if (!DgDx_out.isEmpty()) block_DgDx->LoadBlockValues(*split_DgDx, (*rowIndex)[i]);

    // Assemble multiple steps on this domain into g and dgdp(0) vectors
    if (g_out.get()) g_out->Update(1.0, *split_g, 1.0);

    // HARDWIRED for g is a scalar
    if (!DgDp_out.isEmpty())
       DgDp_out.getMultiVector()->Update(1.0, *((*split_DgDp)(0)), 1.0);

  } // End loop over multiPoint steps on this domain/cluster

  //Copy block vectors into *_out vectors of same size
  if (f_out.get()) f_out->operator=(*block_f);
  if (!DfDp_out.isEmpty()) 
    DfDp_out.getMultiVector()->operator=(*block_DfDp);
  if (!DgDx_out.isEmpty()) 
    DgDx_out.getMultiVector()->operator=(*block_DgDx);

  //Sum together obj fn contributions.
    // HARDWIRED for g is a scalar
  if (numTimeDomains > 1) {
    if (g_out.get()) {
      double g_dist[1];
      double g_sum[1];
      if (globalComm->SubDomainComm().MyPID()==0) g_dist[0] = g_out->operator[](0);
      else g_dist[0]=0;
      globalComm->SumAll(g_dist, g_sum, 1);
      g_out->operator[](0) = g_sum[0];
    }
    if (!DgDp_out.isEmpty()) {
      int  num_p0 =  underlyingME->get_p_map(0)->NumMyElements();
      double* g_dist = new double[num_p0];
      double* g_sum = new double[num_p0];
      for (int i=0; i<num_p0; i++) {
        if (globalComm->SubDomainComm().MyPID()==0)
          g_dist[i] = DgDp_out.getMultiVector()->operator()(0)->operator[](i);
        else g_dist[i]=0;
      }
      globalComm->SumAll(g_dist, g_sum, num_p0);
      for (int i=0; i<num_p0; i++) {
        DgDp_out.getMultiVector()->operator()(0)->operator[](i) = g_sum[i];
      }
      delete g_dist;
      delete g_sum;
    }
  }
}
