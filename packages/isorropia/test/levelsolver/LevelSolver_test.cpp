// includes for the test
#include <iostream>
#include <Epetra_Version.h>
#include <Epetra_SerialComm.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_Time.h>
#include <vector>

#include "LevelSolver.h"
#include "TestUtils.h"
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

#include <Kokkos_SerialNode.hpp>
#ifdef KOKKOS_HAVE_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef KOKKOS_HAVE_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif

#define USE_ISORROPIA

#ifdef USE_ISORROPIA
  #define FNARGS "[matrix filename]"
#else
  #define FNARGS "[matrix filename] [permutation filename]"
#endif

#define PRINT_AND_EXIT() \
{ \
    cout << "Usage: " << argv[0] << " " << FNARGS << " [thread args] [numTrials]" << endl \
                      << "where [thread args] takes one of the following forms:" << endl \
                      << " b:e+i    e.g., 1:5+1 runs 1,2,3,4,5 threads" << endl \
                      << " b:e*m    e.g., 1:2:16 runs 1,2,4,8,16 threads" << endl \
                      << " n[,n]*   e.g., 1,3,9,12 runs 1,3,9,12 threads" << endl; \
    return -1; \
}

using std::endl;
using std::cout;
using std::string;
using std::setw;
using std::setprecision;
using std::fixed;
using std::scientific;
using std::vector;
using std::ifstream;
using std::ofstream;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main() routines for the test
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);

  // cout << Epetra_Version() << endl << endl;
  string mfn;
  vector<int> threads;
  int numTrials = 30;

#ifndef USE_ISORROPIA
  string pfn;
  if (argc > 2) {
    mfn = argv[1];
    pfn = argv[2];
    if (argc > 3) if ( parseNumThreads(argv[3],threads) != 0 ) PRINT_AND_EXIT();
    if (argc > 4) if ( sscanf(argv[4],"%d",&numTrials)  != 1 ) PRINT_AND_EXIT();
  }
  else PRINT_AND_EXIT();
#else
  if (argc > 2) {
    mfn = argv[1];
    if (argc > 2) if ( parseNumThreads(argv[2],threads) != 0 ) PRINT_AND_EXIT();
    if (argc > 3) if ( sscanf(argv[3],"%d",&numTrials)  != 1 ) PRINT_AND_EXIT();
  }
  else PRINT_AND_EXIT();
#endif
  Epetra_CrsMatrix *L;
  int ierr = EpetraExt::MatrixMarketFileToCrsMatrix(mfn.c_str(),Comm,L,false,true); // transpose=false, verbose=true
  if (ierr) {
    cout << "EpetraExt::MatrixMarketFIleToCrsMatrix returned error " << ierr << "." << endl;
    return -1;
  }

  typedef Kokkos::TBBNode Node;
  Node node(threads[0]); threads.resize(1);
  // typedef Kokkos::SerialNode Node;
  // Node node;
  //typedef Kokkos::TPINode Node;
  //Node node;

  typedef Node::buffer<double>::buffer_t dbuffer;
  typedef Node::buffer<int>::buffer_t ibuffer;
  typedef Kokkos::MultiVector<double,int,Node> MV;
  typedef Kokkos::CrsMatrix<double,int,Node> MAT;

  const int NLRs = L->RowMap().NumMyPoints();
  const Kokkos::size_type NNZ = L->NumMyNonzeros();
  Kokkos::size_type *NNZperRow = new Kokkos::size_type[NLRs];
  double stddev = 0;
  double mean = (double)(NNZ) / (double)(NLRs);
  for (int i=0; i<NLRs; ++i) {
    NNZperRow[i] = L->NumMyEntries(i);
    double tmp = (NNZperRow[i] - mean);
    stddev += (tmp*tmp);
  }
  stddev = sqrt(stddev / NLRs);
  MAT Lpar(node);
  Lpar.initializeProfile(NLRs,NNZperRow);
  delete [] NNZperRow; NNZperRow = 0;
  // put matrix data into Kokkos::CrsMatrix
  for (int r=0; r<NLRs; ++r) {
    double *vals;
    int    *inds;
    int     numentries;
    L->ExtractMyRowView(r,numentries,vals,inds);
    Lpar.insertEntries(r,numentries,inds,vals);
  }
  Kokkos::DefaultSparseMultiply<MAT,MV> DSMV(node);
  DSMV.initializeStructure(Lpar,true);
  DSMV.initializeValues(Lpar,true);

  cout << endl << "*** Matrix statistics: " << mfn << endl;
  cout << "Number of rows: " << NLRs << endl;
  cout << "Number of non-zeros: " << NNZ << endl;
  cout << "Mean number of non-zeros per row: " << fixed << setprecision(1) << mean << endl;
  cout << "Std dev number of non-zeros per row: " << fixed << setprecision(2) << stddev << endl;

  Epetra_LevelSolver<Node> LS(L->RowMap(),node);

#ifndef USE_ISORROPIA
  {
    vector<int> pvec(L->NumMyRows());
    int NumLevels;
    vector<int> lsizes;
    try {
      ifstream Pfn;
      Pfn.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit ); 
      Pfn.open(pfn.c_str());
      for (int i=0; i<L->NumMyRows(); ++i) {
        Pfn >> pvec[i];
      }
      Pfn >> NumLevels;
      lsizes.resize(NumLevels);
      for (int i=0; i<NumLevels; ++i) {
        Pfn >> lsizes[i];
      }
    }
    catch (ifstream::failure e) {
      cout << "Exception opening/reading file " << pfn << endl;
      return -1;
    }
    LS.SetLevelInfo(NumLevels,&lsizes[0],&pvec[0]); 
  }
#else
  {
    double time;
    timer.ResetStartTime();
    ierr = LS.Analyze(L->Graph()); 
    time = timer.ElapsedTime();
    cout << "\nLevelSolver::Analyze() time: " << time << endl;
  }
  if (ierr) {
    cout << "LevelSolver::Analyze returned an error" << endl;
    return -1;
  }
#endif
  {
    double time;
    timer.ResetStartTime();
    LS.Setup(*L);
    time = timer.ElapsedTime();
    cout << "LevelSolver::Setup() time: " << time << endl;
  }
  cout << "\n*** LevelSolver statistics" << endl;
  LS.Print(cout,1);
  {
    ofstream fout("levelinfo.dat");
    LS.Print(fout,2);
  }

  // test LevelSolver
  Epetra_Vector x(L->RowMap(),false); 
  Epetra_Vector Lx(x);

  cout << "\n*** Performing verification" << endl;
  for (int t=0; t<4; ++t) {
    x.SetSeed(static_cast<unsigned int>(100000*timer.WallTime()));
    x.Random();
    Lx.Random();
    // cout << "Actual Solution: " << endl; x.Print(cout);
    if ((t & 2) == 0) {
      cout << "Applying L to x using CrsMatrix\n";
      L->Apply(x,Lx);
    }
    else {
      cout << "Applying L to x using LevelSolver\n";
      Lx = x;
      LS.ApplyInverse(Lx,Lx);
    }
    // cout << "RHS: " << endl; Lx.Print(cout);
    if ((t & 1) == 0) {
      cout << "Solving L*x using CrsMatrix\n";
      L->Solve(false,false,false,Lx,Lx);
    }
    else {
      cout << "Solving L*x using LevelSolver\n";
      LS.Apply(Lx,Lx);
    }
    double errnrm, xnrm;
    Lx.Update(-1.0,x,1.0);
    Lx.Norm2(&errnrm);
    x.Norm2(&xnrm);
    cout << "||x - inv(L)*(L*x)||/||x||: " << setprecision(2) << scientific << errnrm/xnrm << "\n\n";
  }

  // for timings, neglect diagonal of LS because the scaling routines are not parallelized
  // LS.setUnitDiag(true);
  LS.setIgnorePerm(true);

  cout << "*** Timings over " << numTrials << " trials" << endl;
  for (vector<int>::iterator nt=threads.begin(); nt != threads.end(); ++nt) {
    double time;
    // node.init(*nt);
    // LevelSolver inverse
    timer.ResetStartTime();
    for (int t=0; t<numTrials; ++t) {
      LS.Apply(Lx,Lx);
    }
    time = timer.ElapsedTime();
    cout << setw(20) << "LevelSolver  solve, " << setw(2) << *nt << " threads, "
         << setprecision(2) << scientific << setw(11) << time           << " total, " 
         << setprecision(2) << scientific << setw(11) << time/numTrials << " average" 
         << endl;
    // LevelSolver forward
    timer.ResetStartTime();
    for (int t=0; t<numTrials; ++t) {
      LS.ApplyInverse(Lx,Lx);
    }
    time = timer.ElapsedTime();
    cout << setw(20) << "LevelSolver  apply, " << setw(2) << *nt << " threads, "
         << setprecision(2) << scientific << setw(11) << time           << " total, " 
         << setprecision(2) << scientific << setw(11) << time/numTrials << " average" 
         << endl;
    // Kokkos::CrsMatrix multiply
    MV x(node), Lx(node);
    dbuffer vecbuf = node.allocBuffer<double>(2*NLRs);
    x.initializeValues(NLRs,1,vecbuf,NLRs);
    Lx.initializeValues(NLRs,1,vecbuf+NLRs,NLRs);
    timer.ResetStartTime();
    for (int t=0; t<numTrials; ++t) {
      DSMV.Apply(false,1.0,x,0.0,Lx);
    }
    time = timer.ElapsedTime();
    node.freeBuffer<double>(vecbuf);
    cout << setw(20) << "K::CrsMatrix apply, " << setw(2) << *nt << " threads, "
         << setprecision(2) << scientific << setw(11) << time            << " total, " 
         << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
         << endl;
  }
  {
    double time;
    // Epetra inverse
    timer.ResetStartTime();
    for (int t=0; t<numTrials; ++t) {
      L->Solve(false,false,false,Lx,Lx);
    }
    time = timer.ElapsedTime();
    cout << setw(20) << "Epetra_CrsMatrix solve,         "
         << setprecision(2) << scientific << setw(11) << time            << " total, " 
         << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
         << endl;
    // Epetra multiply
    timer.ResetStartTime();
    for (int t=0; t<numTrials; ++t) {
      L->Apply(x,Lx);
    }
    time = timer.ElapsedTime();
    cout << setw(20) << "Epetra_CrsMatrix apply,         "
         << setprecision(2) << scientific << setw(11) << time            << " total, " 
         << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
         << endl;
  }

  // delete manually allocated matrices
  delete L;
  return 0;
}
