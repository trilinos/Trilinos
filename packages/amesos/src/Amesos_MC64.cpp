#include "Amesos_ConfigDefs.h"
#include "Amesos_MC64.h"
#include "Epetra_Comm.h"
#include "Epetra_RowMatrix.h"
#include <map>
using namespace std;

extern "C" void F77_FUNC(mc64id, MC64ID)(int*);
extern "C" void F77_FUNC(mc64ad, MC64AD)(int*, int*, int*, int*, int*, 
                                         double*, int*, int*, int*, int*, 
                                         int*, double*, int*, int*);

// ===========================================================================
Amesos_MC64::Amesos_MC64(const Epetra_RowMatrix& A, int JOB,
                         const bool StoreTranspose, const bool analyze) :
  A_(A)
{
  if (A_.Comm().NumProc() != 1)
  {
    cerr << "Class Amesos_MC64 can be used with one processor only!" << endl;
    exit(EXIT_FAILURE);
  }
  F77_FUNC(mc64id, MC64ID)(ICNTL_);

  Compute(JOB, StoreTranspose, analyze);
}

// ===========================================================================
int Amesos_MC64::Compute(int JOB, const bool StoreTranspose, const bool analyze)
{
  // convert A_ into column-based format

  int MaxNumEntries = A_.MaxNumEntries();
  int N = A_.NumMyRows();
  int NE = A_.NumMyNonzeros();
  vector<int> IP;
  vector<int> IRN;
  vector<double> VAL;

  vector<int> Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  if (StoreTranspose)
  {
    // cheapest way, just store the transpose of A and not A. This is because
    // we have easy access to rows and not to columns.
    IP.resize(N + 1); IP[0] = 1;
    IRN.resize(NE);
    VAL.resize(NE);
    int count = 0;
    
    for (int i = 0 ; i < N ; ++i)
    {
      int NumEntries = 0;

      A_.ExtractMyRowCopy(i, MaxNumEntries, NumEntries, &Values[0], 
                          &Indices[0]);

      IP[i + 1] = IP[i] + NumEntries;

      for (int j = 0 ; j < NumEntries ; ++j)
      {
        IRN[count] = Indices[j] + 1;
        VAL[count] = Values[j];
        ++count;
      }
    }
    assert(count == NE);
  }
  else
  {
    // this stores the matrix and not its transpose, but it is more memory 
    // demading. The ifdef'd out part is simple and fast, but very memory
    // demanding.

#if 0
    IRN.resize(N * MaxNumEntries);
    VAL.resize(N * MaxNumEntries);

    vector<int> count(N);
    for (int i = 0 ; i < N ; ++i) count[i] = 0;

    for (int i = 0 ; i < N ; ++i)
    {
      int NumEntries = 0;

      A_.ExtractMyRowCopy(i, MaxNumEntries, NumEntries, &Values[0], 
                          &Indices[0]);

      for (int j = 0 ; j < NumEntries ; ++j)
      {
        int col = Indices[j];
        IRN[col * MaxNumEntries + count[col]] = i + 1;
        VAL[col * MaxNumEntries + count[col]] = Values[j];
        ++count[col];
      }
    }

    // now compact storage
    int k = 0;
    for (int col = 0 ; col < N ; ++col)
    {
      for (int row = 0 ; row < count[col] ; ++row)
      {
        IRN[k] = IRN[col * MaxNumEntries + row];
        VAL[k] = VAL[col * MaxNumEntries + row];
        ++k;
      }
    }
    assert (k == NE);

    IRN.resize(k);
    VAL.resize(k);

    IP.resize(N + 1);
    IP[0] = 1;

    for (int col = 0 ; col < N ; ++col)
      IP[col + 1] = IP[col] + count[col];
#else
    vector<vector<int> > cols(N);
    vector<vector<double> > vals(N);

    for (int i = 1 ; i <= N ; ++i)
    {
      int NumEntries = 0;

      A_.ExtractMyRowCopy(i - 1, MaxNumEntries, NumEntries, &Values[0], 
                          &Indices[0]);

      for (int j = 0 ; j < NumEntries ; ++j)
      {
        cols[Indices[j]].push_back(i);
        vals[Indices[j]].push_back(Values[j]);
      }
    }

    IP.resize(N + 1); IP[0] = 1;
    IRN.resize(NE);
    VAL.resize(NE);
    int count = 0;

    for (int i = 0 ; i < N ; ++i)
    {
      IP[i + 1] = IP[i] + cols[i].size();

      for (int j = 0 ; j < cols[i].size() ; ++j)
      {
        IRN[count] = cols[i][j];
        VAL[count] = vals[i][j];
        ++count;
      }
    }
#endif
  }

  int NUM;
  CPERM_.resize(N);
  int LIW = 10 * N + NE;
  vector<int> IW(LIW);
  int LDW = 3 * N + NE;
  DW_.resize(LDW);

  JOB = 5;
  F77_FUNC(mc64ad, MC64aD)(&JOB, &N, &NE, &IP[0], &IRN[0], &VAL[0], 
                           &NUM, &CPERM_[0], &LIW, &IW[0], &LDW, 
                           &DW_[0], ICNTL_, INFO_);

  if (analyze)
  {
    map<double, int> table;
    for (int col = 0 ; col < N ; ++col)
    {
      for (int j = IP[col] ; j < IP[col + 1] ; ++j)
      {
        int row = IRN[j - 1] - 1;
        int new_col = CPERM_[col];
        double new_val = VAL[j - 1] * exp(DW_[row] + DW_[col + N]);
        if (new_val < 0.0) new_val = -new_val;
        if (new_val > 0.1) table[0.1]++;
        else if (new_val > 0.01) table[0.01]++;
        else if (new_val > 0.001) table[0.001]++;
        else if (new_val > 0.0001) table[0.0001]++;
        else if (new_val > 0.00001) table[0.00001]++;
        else if (new_val > 0.000001) table[0.000001]++;
        else table[0.0]++;
      }
    }

    cout << "# elements (total)    = " << NE << endl;
    cout << "# elements > 0.1      = " << table[0.1] << endl;
    cout << "# elements > 0.01     = " << table[0.01] << endl;
    cout << "# elements > 0.001    = " << table[0.001] << endl;
    cout << "# elements > 0.0001   = " << table[0.0001] << endl;
    cout << "# elements > 0.00001  = " << table[0.00001] << endl;
    cout << "# elements > 0.000001 = " << table[0.000001] << endl;
    cout << "# elements <=0.000001 = " << table[0.0] << endl;
  }

  AMESOS_RETURN(INFO_[0]);
}

// ===========================================================================
double* Amesos_MC64::GetColScaling()
{
  return((double*)&DW_[0 + A_.NumMyRows()]);
}
