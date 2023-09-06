// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER
//

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziTpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziSolverUtils.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_LAPACK.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Tpetra::CrsMatrix;
using Tpetra::Map;

// Function to print a vector for debugging purposes
template <class ST>
void print(std::vector<ST> const &a) {
  std::cout << "The vector elements are : ";

  for(int i=0; i < a.size(); i++)
    std::cout << a.at(i) << ' ';
}

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  typedef double                              ST;
  typedef int                                 OT;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Tpetra::MultiVector<ST>             MV;
  //typedef MV::global_ordinal_type             GO;
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  typedef Anasazi::SolverUtils<ST, MV, OP> msutils;
  const ST ONE  = SCT::one();

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  const int MyPID = comm->getRank ();
  //const int NumImages = comm->getSize ();

  bool testFailed = true;
  //bool verbose = true;
  std::string ortho("ICGS");
  std::string filename("simple.mtx");
  int nev = 4;
  int nsteps = 50;
  int blockSize = 5;
  MT tol = 1.0e-1;
  int res_freq = nsteps+1;
  int ortho_freq = 1;
  int i; 	/* Loop variable */

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("ortho",&ortho,"Orthogonalization method (DGKS, ICGS, or SVQB)");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("nsteps",&nsteps,"Number of times to apply A. ('q' parameter)");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  cmdp.setOption("filename",&filename,"Filename for Matrix Market test matrix.");
  cmdp.setOption("res_freq",&res_freq,"How many iterations between the computation of residuals (Orthogonalization is also done).");
  cmdp.setOption("ortho_freq",&ortho_freq,"How many iterations between the orthogonalization of the basis (Not including checking residuals).");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if ( blockSize < nev ){
    if(MyPID == 0) std::cout << "Block size must be greater than or equal to num evals. Increasing block size to " << nev << std::endl;
    blockSize = nev;
  }
  if (res_freq <= 0) {
    if(MyPID == 0) std::cout << "res_freq must be greater than or equal to 1. Changing res_freq to " << nsteps+1 << std::endl;
    res_freq = nsteps+1;
  }
  if (ortho_freq <= 0) {
    if(MyPID == 0) std::cout << "ortho_freq must be greater than or equal to 1. Changing ortho_freq to " << nsteps+1 << std::endl;
    ortho_freq = nsteps+1;
  }

  //const int ROWS_PER_PROC = 10;
  //int dim = ROWS_PER_PROC * NumImages;
  RCP<CrsMatrix<ST>> A;
  A = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(filename,comm);
  RCP<const Tpetra::Map<> > map = A->getDomainMap();

  // Setup parameters for the eigensolve and residual computation
  RCP<MV> EigenVecs = rcp(new MV(map,blockSize));			/* Store eigenvectors */
  Teuchos::SerialDenseMatrix<OT,ST> H (blockSize, blockSize);		/* Store the Rayleigh-Ritz (RR) problem: H=Q'AQ */
  Teuchos::LAPACK<OT,ST> lapack;					/* Set up Lapack manager */
  Teuchos::SerialDenseMatrix<OT,ST> evects (blockSize, blockSize);	/* Store the right eigenvectors from the RR problem */
  std::vector<MT> evals_real(blockSize);				/* Store the real components of the eigenvalues from RR problem */
  std::vector<MT> evals_imag(blockSize);				/* Store the imag components of the eigenvalues from RR problem */
  std::vector<int> perm(blockSize);					/* Hold the permutation array when sorting the eigenpairs */
  Anasazi::BasicSort<MT> sm;						/* Basic sorting manager */
  int info = -1000;							/* Variable passed into GEEV to check if it returns an error */ 
  ST* vlr = 0; 								/* Placeholder variable that is never accessed */
  const int ldv = 1; 							/* Placeholder variable that is never accessed */
  int lwork = -1;							/* Set to -1 to find the needed size of the workspace */
  std::vector<ST> work(1);						/* Workspace for GEEV */
  std::vector<MT> rwork(2*blockSize);					/* Placeholder variable that is never accessed */
  Anasazi::Eigensolution<ST,MV> sol;					/* Store the solutions to the eigenvalue problem */
  std::vector<Anasazi::Value<ST>> EigenVals(blockSize);			/* Store the eigenvalue solutions */
  RCP<MV> evecs;							/* Store the eigenvectors */
  int numev;								/* Number of eigenpairs found */
  std::vector<MT> normV( blockSize );
  Teuchos::SerialDenseMatrix<int,ST> T (blockSize, blockSize);
  RCP<MV> Avecs;

  // Setting the output format
  std::ostringstream os;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(6);

  // Set up Orthomanager and orthonormalize random vecs. 
  Teuchos::RCP<Anasazi::OrthoManager<ST,MV> > orthoMgr;

  if (ortho == "SVQB") {
      orthoMgr = Teuchos::rcp( new Anasazi::SVQBOrthoManager<ST,MV,OP>());
  } else if (ortho == "DGKS") {
      orthoMgr = Teuchos::rcp( new Anasazi::BasicOrthoManager<ST,MV,OP>());
  } else if (ortho == "ICGS") {
      orthoMgr = Teuchos::rcp( new Anasazi::ICGSOrthoManager<ST,MV,OP>());
  } else {
      TEUCHOS_TEST_FOR_EXCEPTION(ortho!="SVQB"&&ortho!="DGKS"&&ortho!="ICGS",std::logic_error,"Anasazi::RandomSolver Invalid orthogonalization type.");
  }

  // Create initial vectors
  RCP<MV> randVecs = rcp(new MV(map,blockSize));
  randVecs->randomize();
  int rank = orthoMgr->normalize(*randVecs);
  if( rank < blockSize ) std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
  //MVT::MvPrint(*randVecs, std::cout);

  // Applying the subspace iteration
  int num_its = 0;
  while(num_its < nsteps && testFailed) {
    OPT::Apply( *A, *randVecs, *randVecs );
    num_its = num_its + 1;
    if(num_its % ortho_freq == 0 || num_its % res_freq == 0) {
      rank = orthoMgr->normalize(*randVecs);
      if( rank < blockSize ) std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
    }
    if (num_its % res_freq == 0) {
      // Solve RR Problem ---------------------------
      // Building H = Q^TAQ. (RR projection) 
      OPT::Apply( *A, *randVecs, *EigenVecs ); //EigenVecs used for temp storage here. 
      MVT::MvTransMv(ONE, *randVecs, *EigenVecs, H);

      // Querying for the needed size of "work"
      lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
      lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ST>::real (work[0])));
      work.resize( lwork );
  
      // Solve for Harmonic Ritz Values:
      lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
      if(info != 0) std::cout << "Warning!! Anasazi::RandomSolver GEEV solve : info = " << info << std::endl;
      lwork = -1;
	
      // Creating sorting manager to sort the eigenpairs
      sm.sort(evals_real, evals_imag, Teuchos::rcpFromRef(perm), blockSize);
      msutils::permuteVectors(perm, evects);

      // Compute the eigenvalues and eigenvectors from the original eigenproblem
      for( i = 0; i < blockSize; i++){
        EigenVals[i].realpart = evals_real[i];
        EigenVals[i].imagpart = evals_imag[i];
      }
      sol.Evals = EigenVals;

      MVT::MvTimesMatAddMv(ONE,*randVecs,evects,0.0,*EigenVecs);
      sol.Evecs = EigenVecs;
      sol.numVecs = blockSize;

      // Extract evects/evals from solution
      evecs = sol.Evecs;
      numev = sol.numVecs;

      // Verify Evec residuals by hand.
      if (numev > 0) {
        for (i = 0; i < numev; ++i) T(i,i) = sol.Evals[i].realpart;
        Avecs = MVT::Clone( *evecs, numev );
        OPT::Apply( *A, *evecs, *Avecs );
        MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
        MVT::MvNorm( *Avecs, normV );

        testFailed = false;
        for (i = 0; i < numev; i++) {
          if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);	
	  //std::cout << "Iteration " << num_its << " Eval " << sol.Evals[i].realpart << " residual(" << i << "): " << normV[i] << std::endl;
          if (normV[i] > tol ) {
            testFailed = true;
	    break;
	  }
        }
      }
    }
  }

  rank = orthoMgr->normalize(*randVecs);
  if( rank < blockSize ) std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;

  // Building H = Q^TAQ. (RR projection) 
  OPT::Apply( *A, *randVecs, *EigenVecs ); //EigenVecs used for temp storage here. 
  MVT::MvTransMv(ONE, *randVecs, *EigenVecs, H);

  // Solve projected eigenvalue problem.

/* --------------------------------------------------------------------------
 * Parameters for DGEEV:
 *  'N' 		= Don't compute left eigenvectors. 
 *  'V' 		= Compute right eigenvectors. 
 *  blockSize 		= Dimension of H (numEvals)
 *  H.values 		= H matrix (Q'AQ)
 *  H.stride 		= Leading dimension of H (numEvals)
 *  evals_real.data() 	= Array to store evals, real parts
 *  evals_imag.data() 	= Array to store evals, imag parts
 *  vlr 		= Stores left evects, so don't need this
 *  ldv			= Leading dimension of vlr
 *  evects 		= Array to store right evects
 *  evects.stride 	= Leading dimension of evects
 *  work 		= Work array
 *  lwork 		= -1 means to query for array size
 *  rwork 		= Not referenced because ST is not complex
 * -------------------------------------------------------------------------- */

  std::cout << "Starting Harm Ritz val solve (q = " << num_its << ")." << std::endl;
  // Querying for the needed size of "work"
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ST>::real (work[0])));
  work.resize( lwork );
  
  // Solve for Harmonic Ritz Values:
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  if(info != 0) std::cout << "Warning!! Anasazi::RandomSolver GEEV solve : info = " << info << std::endl;
  
  // Creating sorting manager to sort the eigenpairs
  sm.sort(evals_real, evals_imag, Teuchos::rcpFromRef(perm), blockSize);
  msutils::permuteVectors(perm, evects);
 
  // Compute the eigenvalues and eigenvectors from the original eigenproblem
  for( i = 0; i < blockSize; i++){
    EigenVals[i].realpart = evals_real[i];
    EigenVals[i].imagpart = evals_imag[i];
  }

  sol.Evals = EigenVals;

  MVT::MvTimesMatAddMv(ONE,*randVecs,evects,0.0,*EigenVecs);
  sol.Evecs = EigenVecs;
  sol.numVecs = blockSize;

   // Extract evects/evals from solution
  evecs = sol.Evecs;
  numev = sol.numVecs;

  // Computing Rayleigh Quotient: UPDATE - they are identical to what is already computed.
  /* 
  RCP<MV> evecsA = MVT::CloneCopy(*evecs);
  OPT::Apply( *A, *evecs, *evecsA );
  std::vector<ST> RQ(blockSize);
  std::vector<ST> RQ_Scale(blockSize);
  MVT::MvDot(*evecsA, *evecs, RQ); 
  MVT::MvDot(*evecs, *evecs, RQ_Scale); 
  for( i = 0; i < blockSize; i++){
    EigenVals[i] = RQ[i]/RQ_Scale[i]; 
  } 
  */

  std::cout << "Verifying Eigenvector residuals" << std::endl;
  // Verify Evec residuals by hand. 
  if (numev > 0) {
    // Compute the direct residual
    for (i = 0; i < numev; ++i) T(i,i) = sol.Evals[i].realpart;
    Avecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *A, *evecs, *Avecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
    MVT::MvNorm( *Avecs, normV );

    os << "Direct residual norms computed in Tpetra_RandSolve_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
       << "----------------------------------------" << endl;
    testFailed = false;
    for (i=0; i<numev; i++) {
      if ( SCT::magnitude(sol.Evals[i].realpart) != SCT::zero() ) {
        normV[i] = SCT::magnitude(normV[i]/sol.Evals[i].realpart);
      }
      os << std::setw(20) << sol.Evals[i].realpart << std::setw(20) << normV[i] << endl;
      if ( normV[i] > tol ) {
        testFailed = true;
      }
    }
    if (MyPID==0) {
      cout << endl << os.str() << endl;
    }
  }

  if (testFailed) {
    if (MyPID==0) {
      cout << "End Result: TEST FAILED" << endl;
    }
    return -1;
  }
  //
  // Default return value
  //
  if (MyPID==0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}
