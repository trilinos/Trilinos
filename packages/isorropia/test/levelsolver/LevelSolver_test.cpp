// includes for the test
#include <iostream>
#include <Epetra_Version.h>
#include <Epetra_SerialComm.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_Time.h>
#include <vector>


#include "LevelSolver.h"
#include "TestUtils.h"


//#include <Tpetra_MapDecl.hpp>
//#include <Tpetra_ConfigDefs.hpp>

#include <Teuchos_DefaultSerialComm.hpp>


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

//#define USE_ISORROPIA

#ifdef USE_ISORROPIA
  #define FNARGS "[matrix filename]"
#else
  #define FNARGS "[matrix filename] [permutation filename]"
#endif

#define PRINT_AND_EXIT() \
{ \
    cout << "Usage: " << argv[0] << " " << FNARGS << " [numThreads] [numTrials]" << endl; \
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

template<class NodeT>
void verify(const Epetra_CrsMatrix *L, const Epetra_LevelSolver<NodeT> &LS);



template <class nodeT>
void copyEpetravToTpetrav(const Epetra_Vector &evector,
                          Tpetra::Vector<double,int,int,nodeT> &tvector);

template <class nodeT>
void copyTpetravToEpetrav(const Tpetra::Vector<double,int,int,nodeT> &tvector,
                          const Epetra_Vector &evector);


template<class NodeT>
void timeLevelSolver(const Epetra_CrsMatrix *L, Epetra_LevelSolver<NodeT> &LS, int ntrials);

void timeOrigSolver(const Epetra_CrsMatrix *L, int ntrials);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main() routines for the test
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);

  // cout << Epetra_Version() << endl << endl;
  string mfn;
  int numTrials = 30;
  int numThreads = 1;

#ifndef USE_ISORROPIA
  string pfn;
  if (argc > 2) {
    mfn = argv[1];
    pfn = argv[2];
    if (argc > 3) if ( sscanf(argv[3],"%d",&numThreads) != 1 ) PRINT_AND_EXIT();
    if (argc > 4) if ( sscanf(argv[4],"%d",&numTrials)  != 1 ) PRINT_AND_EXIT();
  }
  else PRINT_AND_EXIT();
#else
  if (argc > 2) {
    mfn = argv[1];
    if (argc > 2) if ( sscanf(argv[2],"%d",&numThreads) != 1 ) PRINT_AND_EXIT();
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

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);
  typedef Kokkos::TPINode Node;

  Teuchos::RCP<Node> node = Teuchos::rcp(new Node(pl));

  //typedef Node::buffer<double>::buffer_t dbuffer;
  //Teuchos::ArrayRCP<double> dbuffer;

  typedef Kokkos::MultiVector<double,Node> MV;
  typedef Kokkos::CrsMatrix<double,Node> MAT;

  /////////////////////////////////////////////////////////////
  const int NLRs = L->RowMap().NumMyPoints();
  size_t NNZ = L->NumMyNonzeros(); 
  size_t *NNZperRow = new size_t[NLRs];
  double stddev = 0;
  double mean = (double)(NNZ) / (double)(NLRs);
  for (int i=0; i<NLRs; ++i) 
  {
     NNZperRow[i] = L->NumMyEntries(i);
     double tmp = (NNZperRow[i] - mean);
     stddev += (tmp*tmp);
  }
  stddev = sqrt(stddev / NLRs);

  delete [] NNZperRow; NNZperRow = 0;

  cout << endl << "*** Matrix statistics: " << mfn << endl;
  cout << "Number of rows: " << NLRs << endl;
  cout << "Number of non-zeros: " << NNZ << endl;
  cout << "Mean number of non-zeros per row: " << fixed << setprecision(1) << mean << endl;
  cout << "Std dev number of non-zeros per row: " << fixed << setprecision(2) << stddev << endl;
  /////////////////////////////////////////////////////////////

  Epetra_LevelSolver<Node> LS(L->RowMap(),node);

#define USE_ISORROPIA
#ifndef USE_ISORROPIA
   vector<int> pvec(L->NumMyRows());
   int NumLevels;
   vector<int> lsizes;
   try 
   {
      ifstream Pfn;
      Pfn.exceptions(ifstream::eofbit | ifstream::failbit | ifstream::badbit ); 
      Pfn.open(pfn.c_str());
      for (int i=0; i<L->NumMyRows(); ++i) 
      {
         Pfn >> pvec[i];
      }
      Pfn >> NumLevels;
      lsizes.resize(NumLevels);
      for (int i=0; i<NumLevels; ++i) 
      {
         Pfn >> lsizes[i];
      }
   }
   catch (ifstream::failure e) 
   {
       cout << "Exception opening/reading file " << pfn << endl;
       return -1;
   }
   LS.SetLevelInfo(NumLevels,&lsizes[0],&pvec[0]); 
#else
   {
      double time;
      timer.ResetStartTime();
      ierr = LS.Analyze(L->Graph()); 
      time = timer.ElapsedTime();
      cout << "\nLevelSolver::Analyze() time: " << time << endl;
   }
   if (ierr) 
   {
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
      cout << "\n*** LevelSolver statistics" << endl;
   }
   LS.Print(cout,1);
   {
     ofstream fout("levelinfo.dat");
     LS.Print(fout,2);
   }

   // Verify that the level solver correctly solves the system
   verify(L,LS);

   // Time the level solver 
//   timeLevelSolver(L,LS,numTrials);

//   timeOrigSolver(L,numTrials);


  // delete manually allocated matrices
  delete L;
  return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class nodeT>
void verify(const Epetra_CrsMatrix *L, const Epetra_LevelSolver<nodeT> &LS)
{
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);

//  tmapRCP->describe(*Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false)) , Teuchos::VERB_EXTREME  );

  ///////////////////////////////////////////////////////////////////
  // test LevelSolver
  ///////////////////////////////////////////////////////////////////
  Epetra_Vector x_e(L->RowMap(),false); 
  Epetra_Vector x2_e(x_e); 
  Epetra_Vector Lx_e(x_e);

  typedef Tpetra::Vector<double,int,int,nodeT> TV;

  TV x_t(LS.getTpetraMap(),false); 
  TV x2_t(x_t);
  TV Lx_t(x_t);

  cout << "\n*** Performing verification" << endl;

  double errnrm, xnrm;
  //////////////////////////////////////////////
  // Test 1: Sanity Check
  //////////////////////////////////////////////
  cout << "Verification test 1" << std::endl;
  x_e.SetSeed(static_cast<unsigned int>(100000*timer.WallTime()));
  x_e.Random();
  Lx_e.SetSeed(static_cast<unsigned int>(100000*timer.WallTime()));
  Lx_e.Random();

  cout << "Applying L to x using CrsMatrix\n";
  L->Apply(x_e,Lx_e); // Lx = L * x

  cout << "Solving L*x using CrsMatrix\n";
  L->Solve(false,false,false,Lx_e,x2_e);      // Lx = L^-1 Lx = x
          
  x2_e.Update(-1.0,x_e,1.0);
  x2_e.Norm2(&errnrm);
  x_e.Norm2(&xnrm);
  cout << "||x - inv(L)*(L*x)||/||x||: " << setprecision(2) << scientific << errnrm/xnrm << "\n\n";
  //////////////////////////////////////////////

  //////////////////////////////////////////////
  // Test 2: Sanity Check 2
  //////////////////////////////////////////////
//   cout << "Verification test 2" << std::endl;
//   x_t.randomize();
//   Lx_t.randomize();

//   cout << "Applying L to x using LevelSolver\n";
//   LS.ApplyInverse(x_t,Lx_t); // Lx = L * x

//   cout << "Solving L*x using LevelSolver\n"; // 
//   LS.Apply(Lx_t,x2_t);   // Lx = L^-1 Lx = x

//   x2_t.Update(-1.0,x_t,1.0);
//   x2_t.Norm2(&errnrm);
//   x_t.Norm2(&xnrm);
//   cout << "||x - inv(L)*(L*x)||/||x||: " << setprecision(2) << scientific << errnrm/xnrm << "\n\n";
  //////////////////////////////////////////////

  //////////////////////////////////////////////
  // Test 3: 
  //////////////////////////////////////////////
  cout << "Verification test 3" << std::endl;
  x_e.Random();
  Lx_e.Random();

  cout << "Applying L to x using CrsMatrix\n";
  L->Apply(x_e,Lx_e); // Lx = L * x

  copyEpetravToTpetrav(Lx_e,Lx_t);
  copyEpetravToTpetrav(x_e,x_t);

  cout << "Solving L*x using LevelSolver\n"; // 
  LS.Apply(Lx_t,x2_t);   // Lx = L^-1 Lx = x
  x2_t.update(-1.0,x_t,1.0);

  errnrm = x2_t.norm2();
  xnrm = x_t.norm2();

  cout << "||x - inv(L)*(L*x)||/||x||: " << setprecision(2) << scientific << errnrm/xnrm << "\n\n";
  //////////////////////////////////////////////

  //////////////////////////////////////////////
  // Test 4
  //////////////////////////////////////////////
//   cout << "Verification test 4" << std::endl;
//   x_t.randomize();
//   Lx_t.randomize();

//   cout << "Applying L to x using LevelSolver\n";
//   LS.ApplyInverse(x_t,Lx_t); // Lx = L * x

//   copyTpetravToEpetrav(Lx_t,Lx_e);
//   copyTpetravToEpetrav(x_t,x_e);

//   cout << "Solving L*x using CrsMatrix\n";
//   L->Solve(false,false,false,Lx_e,x2_e);      // Lx = L^-1 Lx = x

//   x2_e.Update(-1.0,x_e,1.0);
//   x2_e.Norm2(&errnrm);
//   x_e.Norm2(&xnrm);
//   cout << "||x - inv(L)*(L*x)||/||x||: " << setprecision(2) << scientific << errnrm/xnrm << "\n\n";
  //////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class nodeT>
void copyEpetravToTpetrav(const Epetra_Vector &evector, 
                          Tpetra::Vector<double,int,int,nodeT> &tvector)
{
  for (int i=0; i< evector.MyLength(); i++)
  {
    tvector.replaceLocalValue(i,evector[i]);
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class nodeT>
void copyTpetravToEpetrav(const Tpetra::Vector<double,int,int,nodeT> &tvector,
		          const Epetra_Vector &evector)
{
  for (int i=0; i< tvector.getLocalLength(); i++)
  {
    evector[i] = tvector.get1dView()[i];
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <class nodeT>
void timeLevelSolver(const Epetra_CrsMatrix *L, Epetra_LevelSolver<nodeT> &LS, int numTrials)
{
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);

  Epetra_Vector x(L->RowMap(),false);
  Epetra_Vector Lx(x);

  // for timings, neglect diagonal of LS because the scaling routines are not parallelized
  // LS.setUnitDiag(true);
  // MMW: Need to parallelize scaling routines
  LS.setIgnorePerm(true);

  cout << "*** Timings over " << numTrials << " trials" << endl;

  //////////////////////////////
  // Level Solver solve
  //////////////////////////////
  double time;
  // node.init(*nt);
  // LevelSolver inverse
  timer.ResetStartTime();
  for (int t=0; t<numTrials; ++t) 
  {
      LS.Apply(Lx,Lx);
  }
  time = timer.ElapsedTime();
  cout << setw(20) << "LevelSolver  solve, " << setw(2) << 0 << " threads, "
       << setprecision(2) << scientific << setw(11) << time           << " total, " 
       << setprecision(2) << scientific << setw(11) << time/numTrials << " average" 
       << endl;
  //////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void timeOrigSolver(const Epetra_CrsMatrix *L, int numTrials)
{
  Epetra_SerialComm Comm;
  Epetra_Time timer(Comm);

  Epetra_Vector x(L->RowMap(),false);
  Epetra_Vector Lx(x);

  // for timings, neglect diagonal of LS because the scaling routines are not parallelized

  cout << "*** Timings over " << numTrials << " trials" << endl;


  //////////////////////////////
  // Original Epetra solver
  //////////////////////////////
  double time;
  // Epetra inverse
  timer.ResetStartTime();
  for (int t=0; t<numTrials; ++t) 
  {
    L->Solve(false,false,false,Lx,Lx);
  }
  time = timer.ElapsedTime();
  cout << setw(20) << "Epetra_CrsMatrix solve,         "
       << setprecision(2) << scientific << setw(11) << time            << " total, " 
       << setprecision(2) << scientific << setw(11) << time /numTrials << " average" 
       << endl;
  //////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////
