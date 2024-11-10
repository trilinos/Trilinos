// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
// This driver reads a matrix from a file, uses the corresponding map to
// create a set of random multivectors, and then orthogonalizes the
// multivectors. The orthogonalization function at the end of this file
// can be included in other codes to orthogonalize an arbitrary
// Tpetra::MultiVector.
//
//
#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_StackedTimer.hpp"

#include "BelosTpetraAdapter.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosTpetraOperator.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "BelosOrthoManagerFactory.hpp"
#include "BelosOrthoManager.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;

//Function Forward Declaration:
template <class ScalarType>
int orthogTpMVecs(Tpetra::MultiVector<ScalarType> & inputVecs, RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> & coeffs, std::string orthogType, int blkSize);

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  {
  typedef double                            ScalarType;
  typedef int                               OT;
  typedef Tpetra::MultiVector<ScalarType>   MV;
  typedef Tpetra::MultiVector<ScalarType>::mag_type   MT;
  typedef Belos::MultiVecTraits<ScalarType,MV>     MVT;
  typedef Teuchos::SerialDenseMatrix<OT,ScalarType> MAT;

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int MyPID = Teuchos::rank(*comm);

  bool verbose = true;
  int blockSize = 4;
  int numVecs = 10;
  int vecLength = 100;
  std::string orthoType("ICGS");
  bool use_stacked_timer = false;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("ortho", &orthoType, "Type of orthogonalization: ICGS, IMGS, DGKS.");
  cmdp.setOption("blkSize",&blockSize,"Number of vectors to orthogonalize at each step. ");
  cmdp.setOption("numVecs",&numVecs,"Total number of vectors for tester to orthogonalize.");
  cmdp.setOption("vecLength",&vecLength,"Length of vectors to be orthogonalized (num elements).");
  cmdp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "Run with or without stacked timer output");

  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  bool proc_verbose = ( verbose && (MyPID==0) ); /* Only print on the zero processor */

  // Create the timer.
  RCP<std::ostream> outputStream = rcp(&std::cout,false);
  RCP<Belos::OutputManager<ScalarType> > printer_ = rcp( new Belos::OutputManager<ScalarType>(Belos::TimingDetails,outputStream) );
  std::string OrthoLabel = "Total Orthog time:";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
    RCP<Teuchos::Time> timerOrtho_ = Teuchos::TimeMonitor::getNewCounter(OrthoLabel);
#endif

  // Set output stream and stacked timer:
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);
  // Set up timers
  Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
  if (use_stacked_timer){
    stacked_timer = rcp(new Teuchos::StackedTimer("Main"));
  }
  Teuchos::TimeMonitor::setStackedTimer(stacked_timer);

  // Create map and random multivec to orthogonalize.
  RCP<const Tpetra::Map<> > map = rcp (new Tpetra::Map<> (vecLength,0,comm));

  RCP<Tpetra::MultiVector<ScalarType>> X1 = rcp( new Tpetra::MultiVector<ScalarType>(map, numVecs) );
  X1->randomize();
  RCP<Tpetra::MultiVector<ScalarType>> XCopy = rcp(new Tpetra::MultiVector<ScalarType>(*X1, Teuchos::Copy)); //Deep copy of X1.

  // Orthogonalize multivec.
  RCP<MAT> coeffMat = rcp(new MAT());
  int rank = -1;
  { //scope guard for timer
  #ifdef BELOS_TEUCHOS_TIME_MONITOR
  Teuchos::TimeMonitor orthotimer(*timerOrtho_);
  #endif
  rank = orthogTpMVecs(*X1, coeffMat, orthoType, blockSize);
  }
  if(proc_verbose){
    std::cout << std::endl << "We have orthogonalized " << numVecs << " vectors. Numerical rank is " << rank << "." << std::endl;
  }

  //Verify Orthogonality:
  RCP<MAT> bigDotAns  = rcp( new MAT(numVecs, numVecs));
  MVT::MvTransMv(1.0, *X1, *X1, *bigDotAns);
  if(proc_verbose){
    std::cout << std::endl << "Printed dot prod matrix for verification: " << std::endl;
    std::cout << "Should be ones on diagonal and zeros elsewhere." << std::endl;
    bigDotAns->print(std::cout);
    std::cout << std::endl << std::endl;
  }

  //Verify coefficients:
  MV coeffs_mv = impl::makeStaticLocalMultiVector (*X1, numVecs, numVecs);
  Tpetra::deep_copy(coeffs_mv, *coeffMat);
  XCopy->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *X1, coeffs_mv, -1.0);
  std::vector<MT> norms(numVecs);
  Teuchos::ArrayView<MT> normView(norms);
  XCopy->norm2(normView);
  if(proc_verbose){
    std::cout << "Here are the QR minus OrigMatrix norms.  Should be zero." << std::endl;
    std::cout << normView << std::endl << std::endl;
  }

  if(proc_verbose){
    //Print final timing details:
    Teuchos::TimeMonitor::summarize( printer_->stream(Belos::TimingDetails) );

    if (use_stacked_timer) {
      stacked_timer->stop("Main");
      Teuchos::StackedTimer::OutputOptions options;
      options.output_fraction = options.output_histogram = options.output_minmax = true;
      stacked_timer->report(out, comm, options);
    }
  }

  }// End Tpetra Scope Guard
  return 0;
}


// This function orthogonalizes the given MultiVector using
// orthogonalization type 'orthogType' and blocks of size
// blkSize.  It returns the corresponding upper-triangular
// matrix of orthogonalization coefficients.
template <class ScalarType>
int orthogTpMVecs(Tpetra::MultiVector<ScalarType> & inputVecs, RCP<Teuchos::SerialDenseMatrix<int,ScalarType>> & coeffs,  std::string orthogType, int blkSize){
  typedef int                               OT;
  typedef typename Teuchos::SerialDenseMatrix<OT,ScalarType> MAT;
  typedef Tpetra::MultiVector<ScalarType>   MV;
  typedef Tpetra::Operator<ScalarType>             OP;
  int numVecs = inputVecs.getNumVectors();

  //Default OutputManager is std::cout.
  Teuchos::RCP<Belos::OutputManager<ScalarType> > myOutputMgr = Teuchos::rcp( new Belos::OutputManager<ScalarType>() );

  //Check that block size is not bigger than total num vecs.
  if(blkSize > numVecs){
    blkSize = numVecs;
  }
  int numLoops = numVecs/blkSize;
  int remainder = numVecs % blkSize;

  coeffs->shape(numVecs,numVecs);
  RCP<MAT> coeffDot;
  std::vector<RCP<const MV>> pastVecArray; // Stores vectors after orthogonalized
  Teuchos::ArrayView<RCP<const MV>> pastVecArrayView;  // To hold the above

  // Create the orthogonalization manager:
  Belos::OrthoManagerFactory<ScalarType, MV, OP> factory;
  Teuchos::RCP<Teuchos::ParameterList> paramsOrtho;   // can be null
  const Teuchos::RCP<Belos::OrthoManager<ScalarType,MV>> orthoMgr =
                    factory.makeOrthoManager (orthogType, Teuchos::null, myOutputMgr, "Tpetra OrthoMgr", paramsOrtho);

  // Get a view of the first block and normalize: (also orthogonalizes these vecs wrt themselves)
  RCP<MV> vecBlock = inputVecs.subViewNonConst(Teuchos::Range1D(0,blkSize-1));
  RCP<MAT> coeffNorm = Teuchos::rcp( new MAT(Teuchos::View, *coeffs, blkSize, blkSize));
  RCP<MV> pastVecBlock;
  int rank = orthoMgr->normalize(*vecBlock, coeffNorm);
  pastVecArray.push_back(vecBlock);

  // Loop over remaining blocks:
  for(int k=1; k<numLoops; k++){
    pastVecArrayView = arrayViewFromVector(pastVecArray);
    vecBlock = inputVecs.subViewNonConst(Teuchos::Range1D(k*blkSize,k*blkSize + blkSize - 1));
    coeffNorm = Teuchos::rcp( new MAT(Teuchos::View, *coeffs, blkSize, blkSize, k*blkSize, k*blkSize));
    coeffDot = Teuchos::rcp(new MAT(Teuchos::View, *coeffs, k*blkSize, blkSize, 0, k*blkSize));
    Array<RCP<MAT>> dotArray; //Tuechos::Array of Matrices for coeffs from dot products
    dotArray.append(coeffDot);
    rank += orthoMgr->projectAndNormalize(*vecBlock, dotArray, coeffNorm, pastVecArrayView);
    pastVecBlock = inputVecs.subViewNonConst(Teuchos::Range1D(0,k*blkSize + blkSize - 1));
    pastVecArray.front() = pastVecBlock;
  }
  if( remainder > 0){
    pastVecArrayView = arrayViewFromVector(pastVecArray);
    vecBlock = inputVecs.subViewNonConst(Teuchos::Range1D(numVecs-remainder, numVecs-1));
    coeffNorm = Teuchos::rcp( new MAT(Teuchos::View, *coeffs, remainder, remainder, numLoops*blkSize, numLoops*blkSize));
    coeffDot = Teuchos::rcp(new MAT(Teuchos::View, *coeffs, numLoops*blkSize, remainder, 0, numLoops*blkSize));
    Array<RCP<MAT>> dotArray; //Tuechos::Array of Matrices for coeffs from dot products
    dotArray.append(coeffDot);
    rank += orthoMgr->projectAndNormalize(*vecBlock, dotArray, coeffNorm, pastVecArrayView);
  }

  return rank;
}
