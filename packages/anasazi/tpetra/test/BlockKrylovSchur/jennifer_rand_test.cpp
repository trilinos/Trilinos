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
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_LAPACK.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Tpetra::CrsMatrix;
using Tpetra::Map;

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
  typedef MV::global_ordinal_type             GO;
  typedef Tpetra::Operator<ST>                OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  const ST ONE  = SCT::one();

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  const int MyPID = comm->getRank ();
  const int NumImages = comm->getSize ();

  bool testFailed;
  bool verbose = true;
  std::string ortho("ICGS");
  std::string filename("bcsstk12.mtx");
  int nev = 4;
  int nsteps = 3;
  int blockSize = 6;
  MT tol = 1.0e-2;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("ortho",&ortho,"Orthogonalization method (DGKS, ICGS, or SVQB)");
  cmdp.setOption("nev",&nev,"Number of eigenvalues to compute.");
  cmdp.setOption("nsteps",&nsteps,"Number of times to apply A. ('q' parameter)");
  cmdp.setOption("blockSize",&blockSize,"Block size for the algorithm.");
  cmdp.setOption("tol",&tol,"Tolerance for convergence.");
  cmdp.setOption("filename",&filename,"Filename for Matrix Market test matrix.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  if ( blockSize < nev ){
    if(MyPID == 0) std::cout << "Block size must be greater than or equal to num evals. Increasing block size to " << nev << std::endl;
    blockSize = nev;
  }

  const int ROWS_PER_PROC = 10;
  int dim = ROWS_PER_PROC * NumImages;
  RCP<CrsMatrix<ST>> A;
  A = Tpetra::MatrixMarket::Reader<CrsMatrix<ST> >::readSparseFile(filename,comm);
  RCP<const Tpetra::Map<> > map = A->getDomainMap();

  // Create initial vectors
  RCP<MV> randVecs = rcp(new MV(map,blockSize));
  randVecs->randomize();
  //std::cout << "Printing orig rand vec: " << std::endl;
  //MVT::MvPrint(*randVecs, std::cout);

  for( int i = 0; i < nsteps; i++ ){
    OPT::Apply( *A, *randVecs, *randVecs );
  }
  //std::cout << "Printing rand vec after " << nsteps << " applies of A." << std::endl;
  //MVT::MvPrint(*randVecs, std::cout);
  //std::cout << "Print randVecs BEFORE orthog:" << std::endl;
  //MVT::MvPrint(*randVecs, std::cout);

  // Set up Orthomanager and orthonormalize random vecs. 
  Teuchos::RCP<Anasazi::OrthoManager<ST,MV> > orthoMgr;
  if (ortho=="SVQB") {
    orthoMgr = Teuchos::rcp( new Anasazi::SVQBOrthoManager<ST,MV,OP>());
  } else if (ortho=="DGKS") {
    //if (ortho_kappa_ <= 0) {
    //  orthoMgr= Teuchos::rcp( new Anasazi::BasicOrthoManager<ST,MV,OP>();
    //} else {
      orthoMgr = Teuchos::rcp( new Anasazi::BasicOrthoManager<ST,MV,OP>());
    //}   
  } else if (ortho=="ICGS") {
    orthoMgr = Teuchos::rcp( new Anasazi::ICGSOrthoManager<ST,MV,OP>());
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(ortho!="SVQB"&&ortho!="DGKS"&&ortho!="ICGS",std::logic_error,"Anasazi::RandomSolver Invalid orthogonalization type.");
  }

  int rank = orthoMgr->normalize(*randVecs);
  if( rank < blockSize ){
    std::cout << "Warning! Anasazi::RandomSolver Random vectors did not have full rank!" << std::endl;
  }
  //std::cout << "Print randVecs after orthog:" << std::endl;
  //MVT::MvPrint(*randVecs, std::cout);

  //Compute H = Q^TAQ. (RR projection) 
  RCP<MV> EigenVecs = rcp(new MV(map,blockSize));
  Teuchos::SerialDenseMatrix<OT,ST> H (blockSize, blockSize);

  OPT::Apply( *A, *randVecs, *EigenVecs ); //EigenVecs used for temp storage here. 
  MVT::MvTransMv(ONE, *randVecs, *EigenVecs, H);
  //std::cout << "printing H: " << std::endl;
  //H.print(std::cout);

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

  //Find workspace size for DGEEV:

  // 'N' no left eigenvectors. 
  // 'V' to compute right eigenvectors. 
  // blockSize = dimension of H
  // H matrix
  // H.stride = leading dimension of H
  // Array to store evals, real parts
  // Array to store evals, imag parts
  // vlr -> stores left evects, so don't need this
  // lead dim of vlr
  // evects =  array to store right evects
  // evects.stride = lead dim ovf evects
  // work = work array
  // lwork -1 means to query for array size
  // rwork - not referenced because ST is not complex
  //std::cout << "Eigenvalse BEFORE LAPACK, real parts are:" << std::endl;
  //print(evals_real);
  //std::cout << std::endl;
  std::cout << "Starting Harm Ritz val solve." << std::endl;
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ST>::real (work[0])));
  work.resize( lwork );
  // Solve for Harmonic Ritz Values:
  lapack.GEEV('N','V',blockSize,H.values(),H.stride(),evals_real.data(),evals_imag.data(),vlr, ldv, evects.values(), evects.stride(), &work[0], lwork, &rwork[0], &info);
  //std::cout << "Resulting eigenvalues are: " << std::endl;
  //print(evals_real);
  //std::cout << std::endl;
  if(info != 0){
    std::cout << "Warning!! Anasazi::RandomSolver GEEV solve : info = " << info << std::endl;
  }
  std::cout << "Past Harm Ritz val solve." << std::endl;
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
  //TODO: Do evals change much if we compute the Rayleigh Quotient for them? ... But we didn't do that in Matlab, did we... didn't use those.





  // Extract evects/evals from solution
  RCP<MV> evecs = sol.Evecs;
  int numev = sol.numVecs;

  std::cout << "Verifying Eigenvector residuals" << std::endl;
  // Verify Evec residuals by hand. 
  if (numev > 0) {
  //std::cout << "In the numev>0 if" << std::endl;
  // Verify Evec residuals by hand. 
    std::ostringstream os;
    os.setf(std::ios::scientific, std::ios::floatfield);
    os.precision(6);

    // Compute the direct residual
    std::vector<MT> normV( numev );
    Teuchos::SerialDenseMatrix<int,ST> T (numev, numev);
    //std::cout << "About to try to access the evals. Bet it crashes here. Haha." << std::endl;
    for (int i = 0; i < numev; ++i) {
      T(i,i) = sol.Evals[i].realpart;
    }
    //std::cout << "If you see this, it didn't really crash there." << std::endl;
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
