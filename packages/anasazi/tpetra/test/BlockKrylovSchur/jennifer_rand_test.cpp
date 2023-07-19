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

  bool testFailed = false;
  //bool verbose = true;
  std::string ortho("ICGS");
  std::string filename("bcsstk12.mtx");
  int nev = 4;
  int nsteps = 3;
  int blockSize = 6;
  MT tol = 1.0e-2;
  bool conv2tol = false;
  int maxsteps = 5*nev;
  int tolstring = nev;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("ortho",&ortho,"Orthogonalization method (DGKS, ICGS, or SVQB)");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("nsteps",&nsteps,"Number of times to apply A. ('q' parameter)");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  cmdp.setOption("filename",&filename,"Filename for Matrix Market test matrix.");
  cmdp.setOption("conv2tol",&conv2tol,"Run until convergence or until maxsteps is reached. (Default: false)");
  cmdp.setOption("maxsteps",&maxsteps,"If conv2tol=true, set the maximum number of steps we can take (Default: 5*nev)");
  cmdp.setOption("tolstride",&stride,"If conv2tol=true, how often should we check the residuals? (Default: nev)");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if ( blockSize < nev ){
    if(MyPID == 0) std::cout << "Block size must be greater than or equal to num evals. Increasing block size to " << nev << std::endl;
    blockSize = nev;
  }

  //const int ROWS_PER_PROC = 10;
  //int dim = ROWS_PER_PROC * NumImages;
  RCP<CrsMatrix<ST>> A;
  A = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(filename,comm);
  RCP<const Tpetra::Map<> > map = A->getDomainMap();

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
  for( int i = 0; i < nsteps; i++ ) {
    OPT::Apply( *A, *randVecs, *randVecs );
    rank = orthoMgr->normalize(*randVecs);
  }

  if( rank < blockSize ) std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;

  // Compute H = Q^TAQ. (RR projection) 
  RCP<MV> EigenVecs = rcp(new MV(map,blockSize));
  Teuchos::SerialDenseMatrix<OT,ST> H (blockSize, blockSize);

  OPT::Apply( *A, *randVecs, *EigenVecs ); //EigenVecs used for temp storage here. 
  MVT::MvTransMv(ONE, *randVecs, *EigenVecs, H);
  // std::cout << "printing H: " << std::endl;
  // H.print(std::cout);

  // Solve projected eigenvalue problem.
  Teuchos::LAPACK<OT,ST> lapack;
  Teuchos::SerialDenseMatrix<OT,ST> evects (blockSize, blockSize);
  std::vector<MT> evals_real(blockSize);
  std::vector<MT> evals_imag(blockSize);

  int info = -1000; 
  ST* vlr = 0; 
  const int ldv = 1; 

  // Size of workspace and workspace for DGEEV
  int lwork = -1;
  std::vector<ST> work(1);
  std::vector<MT> rwork(2*blockSize);

/* --------------------------------------------------------------------------
 * Parameters for DGEEV:
  // 'N' 		= Don't compute left eigenvectors. 
  // 'V' 		= Compute right eigenvectors. 
  // blockSize 		= Dimension of H (numEvals)
  // H.values 		= H matrix (Q'AQ)
  // H.stride 		= Leading dimension of H (numEvals)
  // evals_real.data() 	= Array to store evals, real parts
  // evals_imag.data() 	= Array to store evals, imag parts
  // vlr 		= Stores left evects, so don't need this
  // ldv		= Leading dimension of vlr
  // evects 		= Array to store right evects
  // evects.stride 	= Leading dimension of evects
  // work 		= Work array
  // lwork 		= -1 means to query for array size
  // rwork 		= Not referenced because ST is not complex
  // */
  std::cout << "Starting Harm Ritz val solve (nsteps = " << nsteps << ")." << std::endl;

  // Querying for the needed size of "work"
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ST>::real (work[0])));
  work.resize( lwork );
  
  // Solve for Harmonic Ritz Values:
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  
  //std::cout << "Resulting eigenvalues are: " << std::endl;
  //print(evals_real);
  //std::cout << std::endl;
  if(info != 0)
    std::cout << "Warning!! Anasazi::RandomSolver GEEV solve : info = " << info << std::endl;
  
  std::cout << "Past Harm Ritz val solve." << std::endl;

  // Creating sorting manager to sort the eigenpairs
  std::vector<int> perm(blockSize);
  Anasazi::BasicSort<MT> sm;
  sm.sort(evals_real, evals_imag, Teuchos::rcpFromRef(perm), blockSize);
  msutils::permuteVectors(perm, evects);
 
  // Compute the eigenvalues and eigenvectors from the original eigenproblem
  Anasazi::Eigensolution<ST,MV> sol;
  std::vector<Anasazi::Value<ST>> EigenVals(blockSize);
  for( int i = 0; i < blockSize; i++){
    EigenVals[i].realpart = evals_real[i];
    EigenVals[i].imagpart = evals_imag[i];
  }

  sol.Evals = EigenVals;

  MVT::MvTimesMatAddMv(ONE,*randVecs,evects,0.0,*EigenVecs);
  sol.Evecs = EigenVecs;
  sol.numVecs = blockSize;

   // Extract evects/evals from solution
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  // Computing Rayleigh Quotient: UPDATE - they are identical to what is already computed.
  /* RCP<MV> evecsA = MVT::CloneCopy(*evecs);
  OPT::Apply( *A, *evecs, *evecsA );
  std::vector<ST> RQ(blockSize);
  std::vector<ST> RQ_Scale(blockSize);
  MVT::MvDot(*evecsA, *evecs, RQ); 
  MVT::MvDot(*evecs, *evecs, RQ_Scale); 
  for( int i = 0; i < blockSize; i++){
    EigenVals[i] = RQ[i]/RQ_Scale[i]; 
  } */

  std::cout << "Verifying Eigenvector residuals" << std::endl;
  // Verify Evec residuals by hand. 
  if (numev > 0) {
    //std::cout << "In the numev>0 if" << std::endl;
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MT> normV( numev );
    Teuchos::SerialDenseMatrix<int,ST> T (numev, numev);
    for (int i = 0; i < numev; ++i) {
      T(i,i) = sol.Evals[i].realpart;
    }
    RCP<MV> Avecs = MVT::Clone( *evecs, numev );

    OPT::Apply( *A, *evecs, *Avecs );

    MVT::MvTimesMatAddMv( -ONE, *evecs, T, ONE, *Avecs );
    MVT::MvNorm( *Avecs, normV );

    os << "Direct residual norms computed in Tpetra_RandSolve_test.exe" << endl
       << std::setw(20) << "Eigenvalue" << std::setw(20) << "Residual  " << endl
       << "----------------------------------------" << endl;
    for (int i=0; i<numev; i++) {
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
