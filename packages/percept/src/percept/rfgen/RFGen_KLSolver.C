#include "RFGen_KLSolver.h"

#include "mpi.h"

// Anasazi basic code
#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziEpetraAdapter.hpp"

// Anasazi solver managers
#include "AnasaziBlockDavidsonSolMgr.hpp"

#include "Epetra_CrsMatrix.h"

namespace RFGen
{

KLSolver::KLSolver(const unsigned spatialDim)
  : 
  Epetra_Operator(),
  mpiComm_(MPI_COMM_WORLD),
  epetraComm_(Epetra_MpiComm(mpiComm_)),
  localNumElem_(0),
  globalNumElem_(0),
  localMaxIntgPts_(0), 
  maxIntgPts_(0), 
  spatialDim_(spatialDim),
  useMatrixFree_(false),
  log_debug_(false)
{}

void 
KLSolver::solve(const int maxNev)
{
  // TODO: template entire class on scalar type
  typedef double ScalarType;
  typedef Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

  appInterface_->computeLocalIntgDataSizes(
    localNumElem_,
    localMaxIntgPts_);

  int myPID;

  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Comm_rank(mpiComm_, &myPID));

  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Allreduce(&localNumElem_, &globalNumElem_, 1, MPI_INT, MPI_SUM, mpiComm_));
  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Allreduce(&localMaxIntgPts_, &maxIntgPts_, 1, MPI_INT, MPI_MAX, mpiComm_));

  // build Epetra objects
  epetra_Map_ = Teuchos::rcp(
    new Epetra_Map(
      globalNumElem_, 
      localNumElem_, 
      0, // IndexBase, 
      epetraComm_));
  
  localElemIDs_.resize(localNumElem_);
  epetra_Map_->MyGlobalElements(&localElemIDs_[0]);

  if (log_debug_)
  {
    epetra_Map_->Print(std::cout);
  }

  // TODO: can we create an "ArrayContainer" class to do this for us??
  std::vector<double> localIntgPtCoords_mem(
    localNumElem_*maxIntgPts_*spatialDim_);
  std::vector<double> localVolumeWeights_mem(
    localNumElem_*maxIntgPts_);

  localIntgPtCoords_.assign(&localIntgPtCoords_mem[0], localNumElem_, maxIntgPts_, spatialDim_);
  localVolumeWeights_.assign(&localVolumeWeights_mem[0], localNumElem_, maxIntgPts_);

  appInterface_->computeLocalIntgData(
    localIntgPtCoords_, 
    localVolumeWeights_);

  epetra_massMat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *epetra_Map_, 1));

  epetra_massMat_->PutScalar(0.0);

  assembleMassMat();

  if (log_debug_)
  {
    epetra_massMat_->Print(std::cout);
  }  

  epetra_massMat_->FillComplete();

  if (useMatrixFree_)
  {
    //assembleMatrixFreeCovariance();
    covarianceOp_ = Teuchos::rcp(this, false);
  }
  else
  {
    assembleFullCovariance();
    covarianceOp_ = epetra_covarianceMat_;
  }

  const int blockSize=pl_.get<int>("Block Size");

  Teuchos::RCP<Epetra_MultiVector> ivec = 
    Teuchos::rcp( new Epetra_MultiVector(*epetra_Map_, blockSize) );
  ivec->PutScalar(0.0);
  
  Teuchos::RCP<Anasazi::BasicEigenproblem<ScalarType,Epetra_MultiVector,Epetra_Operator> > problem =
    Teuchos::rcp( new Anasazi::BasicEigenproblem<ScalarType,Epetra_MultiVector,Epetra_Operator>
		  (covarianceOp_,
		   epetra_massMat_,
		   ivec) );
  
  problem->setHermitian(true);
  
  problem->setNEV(maxNev);
  
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!problem->setProblem(), 
		      "Anasazi Eigenproblem initialized incorrectly");
  
  // TODO: use sublist for solver params
  Anasazi::BlockDavidsonSolMgr<ScalarType,Epetra_MultiVector,Epetra_Operator> MySolverMgr(problem, pl_);
  
  Anasazi::ReturnType returnCode = MySolverMgr.solve();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(returnCode != Anasazi::Converged,
		      "Anasazi Eigensolver failed to converge");
  
  if (0==myPID) {
    std::cout << "Number of iterations = " << MySolverMgr.getNumIters() << std::endl;
  }
  
  // TODO: compute stochastic dimension
  const int stochastic_dim = maxNev;
  TEUCHOS_TEST_FOR_EXCEPT(stochastic_dim > maxNev);

  // write KL solution back to app
  Anasazi::Eigensolution<ScalarType,Epetra_MultiVector> sol = problem->getSolution();

  lambda_mem_.resize(stochastic_dim);
  shards::Array<double,shards::NaturalOrder,Eigen> lambda(&lambda_mem_[0], stochastic_dim);

  for (int i=0; i<stochastic_dim; i++)
  {
    lambda[i] = sol.Evals[i].realpart;
  }
  
  // normalize phi vectors against non-const random vector  
  srand(345608394);

  std::vector<double> coeff(spatialDim_+1);
  for (int d=0; (unsigned)d<spatialDim_+1; d++)
    coeff[d] = (double)rand()/(double)RAND_MAX;

  Teuchos::RCP<Epetra_Vector> test_vec = 
    Teuchos::rcp( new Epetra_Vector(*epetra_Map_) );
  
  for (int i=0; i<localNumElem_; i++)
  {
    double avg_test_vec = 0.0;
    double area = 0.0;

    shards::Array<double,shards::NaturalOrder,Point> JxW1(
      &localVolumeWeights_[i*maxIntgPts_], 
      maxIntgPts_);

    for (int q1=0; q1<maxIntgPts_; q1++)
    {
      double integrand = coeff[0];

      for (int d=0; (unsigned)d<spatialDim_; d++)
	integrand += coeff[d+1]*localIntgPtCoords_(i,q1,d);

      avg_test_vec += integrand * JxW1(q1);
      area += JxW1(q1);
    }

    (*test_vec)[i] = avg_test_vec / area;
  }

  // normalize sign of eigenfunctions: require that sum of components >= 0
  const double tolerance = pl_.get<double>("Convergence Tolerance");
  for (int i=0; i<stochastic_dim; i++)
  {
    double sum = 0;

    Epetra_Vector * phi_vec = (*sol.Evecs)(i);
  
    test_vec->Dot(*phi_vec, &sum);

    if ( (std::fabs(sum) > tolerance) && (sum < 0.0) )
    { 
      double scale = -1.0;    
      phi_vec->Scale(scale);
    }
  }
    
  // generate a non-persisting ArrayView of phi
  // app must copy this data to its native data structures
  const Epetra_MultiVector * phi_vecs = &*(sol.Evecs);
  shards::Array<double,shards::NaturalOrder,Eigen,Cell> phi(phi_vecs->Values(), stochastic_dim, localNumElem_);

  appInterface_->setKLSolution(
    stochastic_dim,
    lambda,
    phi);

  // release local memory
  localIntgPtCoords_mem.clear();
  localVolumeWeights_mem.clear();
}

void
KLSolver::assembleFullCovariance()
{
  epetra_covarianceMat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *epetra_Map_, globalNumElem_));

  epetra_covarianceMat_->PutScalar(0.0);

  // 1) local assembly of covar mat: diagonal square blocks
  {
    int nonlocalNumElem = localNumElem_;
    nonlocalElemIDs_ = localElemIDs_;

    shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> nonlocalIntgPtCoords(
      localIntgPtCoords_.contiguous_data(), 
      nonlocalNumElem, maxIntgPts_, spatialDim_);

    shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights(
      localVolumeWeights_.contiguous_data(), 
      nonlocalNumElem, maxIntgPts_);

    partialAssemble(nonlocalIntgPtCoords, nonlocalVolumeWeights);
  }

  // 2) parallel assembly of covar mat: off diagonal rectangular blocks
  int numProcs, myPID;
  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Comm_size(mpiComm_, &numProcs)); 
  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Comm_rank(mpiComm_, &myPID));

  MPI_Status status;

  std::vector<double> nonlocalIntgPtCoords_mem;
  std::vector<double> nonlocalVolumeWeights_mem;

  // loop over columns, from right to left, from the diagonal
  for (int ncol=1; ncol<numProcs; ncol++)
  {
    const int rec_from_colPID = (myPID + ncol) % numProcs;
    const int send_to_colPID = (myPID + numProcs - ncol) % numProcs;
    
    // send/receive local num elements
    int nonlocalNumElem = 0;

    const int mpi_tag_numelem = 0;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(&localNumElem_, 1, MPI_INT, send_to_colPID, mpi_tag_numelem, 
						&nonlocalNumElem, 1, MPI_INT, rec_from_colPID, mpi_tag_numelem, 
						mpiComm_, &status));

    // s/r global epetra map ids
    nonlocalElemIDs_.resize(nonlocalNumElem);

    const int mpi_tag_elemids = 1;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(&localElemIDs_[0], localNumElem_, MPI_INT, 
						send_to_colPID, mpi_tag_elemids, 
						&nonlocalElemIDs_[0], nonlocalNumElem, MPI_INT, 
						rec_from_colPID, mpi_tag_elemids, 
						mpiComm_, &status));

    // s/r intg pt coords
    const int size_localIntgPtCoords = localNumElem_*maxIntgPts_*spatialDim_;
    const int size_nonlocalIntgPtCoords = nonlocalNumElem*maxIntgPts_*spatialDim_;
    nonlocalIntgPtCoords_mem.resize(size_nonlocalIntgPtCoords);

    shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> nonlocalIntgPtCoords(
      &nonlocalIntgPtCoords_mem[0], 
      nonlocalNumElem, maxIntgPts_, spatialDim_);

    const int mpi_tag_coords = 2;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(localIntgPtCoords_.contiguous_data(), size_localIntgPtCoords, MPI_DOUBLE, 
						send_to_colPID, mpi_tag_coords, 
						nonlocalIntgPtCoords.contiguous_data(), size_nonlocalIntgPtCoords, MPI_DOUBLE, 
						rec_from_colPID, mpi_tag_coords, 
						mpiComm_, &status));

    // s/r volume weights
    const int size_localVolumeWeights = localNumElem_*maxIntgPts_;
    const int size_nonlocalVolumeWeights = nonlocalNumElem*maxIntgPts_;
    nonlocalVolumeWeights_mem.resize(size_nonlocalVolumeWeights);

    shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights(
      &nonlocalVolumeWeights_mem[0], 
      nonlocalNumElem, maxIntgPts_);

    const int mpi_tag_vw = 3;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(localVolumeWeights_.contiguous_data(), size_localVolumeWeights, MPI_DOUBLE, 
						send_to_colPID, mpi_tag_vw, 
						nonlocalVolumeWeights.contiguous_data(), size_nonlocalVolumeWeights, MPI_DOUBLE, 
						rec_from_colPID, mpi_tag_vw, 
						mpiComm_, &status));

    partialAssemble(nonlocalIntgPtCoords, nonlocalVolumeWeights);
  }

  epetra_covarianceMat_->FillComplete();

  if (log_debug_)
  {
    epetra_covarianceMat_->Print(std::cout);
  }  
}

void 
KLSolver::assembleMassMat()
{
  double localMass;

  for (int row=0; row<localNumElem_; row++)
  {
    shards::Array<double,shards::NaturalOrder,Point,Dim> x1(
      &localIntgPtCoords_[row*maxIntgPts_*spatialDim_], 
      maxIntgPts_,
      spatialDim_);

    shards::Array<double,shards::NaturalOrder,Point> JxW1(
      &localVolumeWeights_[row*maxIntgPts_], 
      maxIntgPts_);

    localMass = 0.0;
    for (int q1=0; q1<maxIntgPts_; q1++)
    {
      localMass += JxW1(q1);
    }

    int RowIndex = localElemIDs_[row];
    int * Indices = &localElemIDs_[row];
    epetra_massMat_->InsertGlobalValues(RowIndex, 1, &localMass, Indices);
  }
}

void 
KLSolver::partialAssemble(
  const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,	       
  const shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights)
{
  double localCovar;
  std::vector<double> covar_mem(maxIntgPts_*maxIntgPts_);
  shards::Array<double,shards::NaturalOrder,Point,Point> covar_values(&covar_mem[0], maxIntgPts_, maxIntgPts_);

  const int nonlocalNumElem = nonlocalIntgPtCoords.dimension(0);

  for (int row=0; row<localNumElem_; row++)
  {
    shards::Array<double,shards::NaturalOrder,Point,Dim> x1(
      &localIntgPtCoords_[row*maxIntgPts_*spatialDim_], 
      maxIntgPts_,
      spatialDim_);

    shards::Array<double,shards::NaturalOrder,Point> JxW1(
      &localVolumeWeights_[row*maxIntgPts_], 
      maxIntgPts_);

    int RowIndex = localElemIDs_[row];
    int * Indices = &localElemIDs_[row];

    for (int col=0; col<nonlocalNumElem; col++)
    {
      shards::Array<double,shards::NaturalOrder,Point,Dim> x2(
        &nonlocalIntgPtCoords[col*maxIntgPts_*spatialDim_], 
	maxIntgPts_,
	spatialDim_);

      shards::Array<double,shards::NaturalOrder,Point> JxW2(
        &nonlocalVolumeWeights[col*maxIntgPts_], 
        maxIntgPts_);

      covarFunc_->computeValues(x1, x2, covar_values);

      localCovar = 0.0;
      for (int q1=0; q1<maxIntgPts_; q1++)
      {
	for (int q2=0; q2<maxIntgPts_; q2++)
	{
	  localCovar += covar_values(q1, q2) * JxW1(q1) * JxW2(q2);
	}
      }

      Indices = &nonlocalElemIDs_[col];
      epetra_covarianceMat_->InsertGlobalValues(RowIndex, 1, &localCovar, Indices);
    }
  }  
}

int 
KLSolver::SetUseTranspose(bool UseTranspose) 
{ 
  /* NO tranpose supported!! */ 
  return -1;
}

int 
KLSolver::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Y = A*X
  
  Y.PutScalar(0.0);
  
  const int numVecs = X.NumVectors();
  TEUCHOS_TEST_FOR_EXCEPT(numVecs != Y.NumVectors()); 

  shards::Array<double,shards::NaturalOrder,Eigen,Cell> Y_vec(Y.Values(), numVecs, localNumElem_);
  shards::Array<double,shards::NaturalOrder,Eigen,Cell> X_vec(X.Values(), numVecs, localNumElem_);

  // 1) local matvec of covar mat: diagonal square blocks
  {
    shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> nonlocalIntgPtCoords(
      localIntgPtCoords_.contiguous_data(), 
      localNumElem_, maxIntgPts_, spatialDim_);

    shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights(
      localVolumeWeights_.contiguous_data(), 
      localNumElem_, maxIntgPts_);

    partialMatVec(nonlocalIntgPtCoords, nonlocalVolumeWeights, X_vec, Y_vec);
  }

  // 2) parallel matvec of covar mat: off diagonal rectangular blocks
  int numProcs, myPID;
  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Comm_size(mpiComm_, &numProcs)); 
  TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Comm_rank(mpiComm_, &myPID));

  if (numProcs==1) return 0;

  MPI_Status status;

  std::vector<double> nonlocalIntgPtCoords_mem;
  std::vector<double> nonlocalVolumeWeights_mem;
  std::vector<double> nonlocalX_vec_mem;

  // loop over columns, from right to left, from the diagonal
  for (int ncol=1; ncol<numProcs; ncol++)
  {
    const int rec_from_colPID = (myPID + ncol) % numProcs;
    const int send_to_colPID = (myPID + numProcs - ncol) % numProcs;
    
    // send/receive local num elements
    int localNumElem = localNumElem_; // avoid passing non-const pointer to MPI_Sendrecv below
    int nonlocalNumElem = 0;

    const int mpi_tag_numelem = 0;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(&localNumElem, 1, MPI_INT, send_to_colPID, mpi_tag_numelem, 
						&nonlocalNumElem, 1, MPI_INT, rec_from_colPID, mpi_tag_numelem, 
						mpiComm_, &status));

    // s/r intg pt coords
    const int size_localIntgPtCoords = localNumElem_*maxIntgPts_*spatialDim_;
    const int size_nonlocalIntgPtCoords = nonlocalNumElem*maxIntgPts_*spatialDim_;
    nonlocalIntgPtCoords_mem.resize(size_nonlocalIntgPtCoords);

    shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> nonlocalIntgPtCoords(
      &nonlocalIntgPtCoords_mem[0], 
      nonlocalNumElem, maxIntgPts_, spatialDim_);

    const int mpi_tag_coords = 1;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(localIntgPtCoords_.contiguous_data(), size_localIntgPtCoords, MPI_DOUBLE, 
						send_to_colPID, mpi_tag_coords, 
						nonlocalIntgPtCoords.contiguous_data(), size_nonlocalIntgPtCoords, MPI_DOUBLE, 
						rec_from_colPID, mpi_tag_coords, 
						mpiComm_, &status));

    // s/r volume weights
    const int size_localVolumeWeights = localNumElem_*maxIntgPts_;
    const int size_nonlocalVolumeWeights = nonlocalNumElem*maxIntgPts_;
    nonlocalVolumeWeights_mem.resize(size_nonlocalVolumeWeights);

    shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights(
      &nonlocalVolumeWeights_mem[0], 
      nonlocalNumElem, maxIntgPts_);

    const int mpi_tag_vw = 2;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(localVolumeWeights_.contiguous_data(), size_localVolumeWeights, MPI_DOUBLE, 
						send_to_colPID, mpi_tag_vw, 
						nonlocalVolumeWeights.contiguous_data(), size_nonlocalVolumeWeights, MPI_DOUBLE, 
						rec_from_colPID, mpi_tag_vw, 
						mpiComm_, &status));

    // s/r X_vec
    const int size_localX_vec = localNumElem_*numVecs;
    const int size_nonlocalX_vec = nonlocalNumElem*numVecs;
    nonlocalX_vec_mem.resize(size_nonlocalX_vec);

    shards::Array<double,shards::NaturalOrder,Eigen,Cell> nonlocalX_vec(
      &nonlocalX_vec_mem[0],
      numVecs, nonlocalNumElem);

    const int mpi_tag_xvec = 3;
    TEUCHOS_TEST_FOR_EXCEPT(MPI_SUCCESS != MPI_Sendrecv(X_vec.contiguous_data(), size_localX_vec, MPI_DOUBLE, 
						send_to_colPID, mpi_tag_xvec, 
						nonlocalX_vec.contiguous_data(), size_nonlocalX_vec, MPI_DOUBLE, 
						rec_from_colPID, mpi_tag_xvec, 
						mpiComm_, &status));

    partialMatVec(nonlocalIntgPtCoords, nonlocalVolumeWeights, nonlocalX_vec, Y_vec);
  }

  return 0;
}

int 
KLSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // NO inverse supported!!
  return -1;
}

double 
KLSolver::NormInf() const
{
  //
  return 0.0;
}

const char * 
KLSolver::Label() const
{
  return "KLSolver";
}

bool 
KLSolver::UseTranspose() const
{
  return false;
}

bool 
KLSolver::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
KLSolver::Comm() const
{
  return epetraComm_;
}

const Epetra_Map & 
KLSolver::OperatorDomainMap() const
{
  return *epetra_Map_;
}

const Epetra_Map & 
KLSolver::OperatorRangeMap() const
{
  return *epetra_Map_;
}

void 
KLSolver::partialMatVec(
  const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,
  const shards::Array<double,shards::NaturalOrder,Cell,Point> &nonlocalVolumeWeights,
  const shards::Array<double,shards::NaturalOrder,Eigen,Cell> &X_vec,
  shards::Array<double,shards::NaturalOrder,Eigen,Cell> &Y_vec) const
{
  if (localNumElem_==0) return;

  const int nonlocalNumElem = nonlocalIntgPtCoords.dimension(0);
  if (nonlocalNumElem==0) return;

  double integral;
  std::vector<double> covar_mem(maxIntgPts_*maxIntgPts_);
  shards::Array<double,shards::NaturalOrder,Point,Point> covar_values(&covar_mem[0], maxIntgPts_, maxIntgPts_);

  const int numVecs = X_vec.dimension(0);

  for (int row=0; row<localNumElem_; row++)
  {
    shards::Array<double,shards::NaturalOrder,Point,Dim> x1(
      &localIntgPtCoords_[row*maxIntgPts_*spatialDim_], 
      maxIntgPts_,
      spatialDim_);

    shards::Array<double,shards::NaturalOrder,Point> JxW1(
      &localVolumeWeights_[row*maxIntgPts_], 
      maxIntgPts_);

    for (int col=0; col<nonlocalNumElem; col++)
    {
      shards::Array<double,shards::NaturalOrder,Point,Dim> x2(
        &nonlocalIntgPtCoords[col*maxIntgPts_*spatialDim_], 
	maxIntgPts_,
	spatialDim_);

      shards::Array<double,shards::NaturalOrder,Point> JxW2(
        &nonlocalVolumeWeights[col*maxIntgPts_], 
        maxIntgPts_);

      covarFunc_->computeValues(x1, x2, covar_values);

      integral = 0.0;
      for (int q1=0; q1<maxIntgPts_; q1++)
      {
	for (int q2=0; q2<maxIntgPts_; q2++)
	{
	  integral += covar_values(q1, q2) * JxW1(q1) * JxW2(q2);
	}
      }
     
      for (int i=0; i<numVecs; i++)
	Y_vec(i,row) += integral * X_vec(i,col);
    }
  }  
}

Teuchos::RCP<KLSolver> 
buildKLSolver(
  const Teuchos::RCP<API_KLSolver> &appInterface,
  const Teuchos::RCP<const CovarianceFunction> &covarFunc,
  const Teuchos::ParameterList &pl,
  const bool useMatrixFree
  )
{
  Teuchos::RCP<KLSolver> klsolver = 
    Teuchos::rcp(new KLSolver(appInterface->getSpatialDim()));

  klsolver->setAPI(appInterface);
  klsolver->setCovarianceFunction(covarFunc);
  klsolver->setPL(pl);
  klsolver->setMatrixFree(useMatrixFree);

  return klsolver;
}

}
