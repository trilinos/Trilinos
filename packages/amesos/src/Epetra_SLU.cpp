// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Epetra_SLU.h"

namespace SLU
{
  //FIXME
  typedef int int_t;
  
extern "C" {
#include "dsp_defs.h"
}
}

#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

struct SLUData
{
  SLU::SuperMatrix A, B, X, L, U;
  SLU::factor_param_t iparam;
  SLU::mem_usage_t mem_usage;
};

Epetra_SLU::Epetra_SLU( Epetra_LinearProblem * Problem,
                        int fill_fac,
                        int panel_size,
                        int relax )
: count_(0)
{
  using namespace SLU;

  data_ = new SLUData();

  data_->iparam.panel_size = panel_size;
  data_->iparam.relax      = relax;
  data_->iparam.diag_pivot_thresh = -1;
  data_->iparam.drop_tol = -1;

  A_ = dynamic_cast<Epetra_CrsMatrix *> (Problem->GetOperator());
  X_ = Problem->GetLHS();
  B_ = Problem->GetRHS();

  AG_ = new Epetra_CrsGraph( A_->Graph() );

  int NumIndices;
  int NumMyCols = AG_->NumMyCols();
  int NumMyEqs = A_->NumMyRows();
  int GlobalMaxNumIndices = AG_->GlobalMaxNumIndices();
  int * xIndices;

  TransNumNZ_ = new int[NumMyCols];
  for( int i = 0; i < NumMyCols; i++ )
    TransNumNZ_[i] = 0;
  for( int i = 0; i < NumMyEqs; i++ )
  {
    assert( AG_->ExtractMyRowView( i, NumIndices, xIndices ) == 0 );
    for( int j = 0; j < NumIndices; j++ )
      ++TransNumNZ_[ xIndices[j] ];
  }
  TotNumNZ_ = 0;
  for( int i = 0; i < NumMyCols; i++ )
   TotNumNZ_ += TransNumNZ_[i];

  TransIndices_ = new int*[ NumMyCols ];
  TransValues_ = new double*[ NumMyCols ];

  RowIndices_ = new int[TotNumNZ_];
  RowValues_ = new double[TotNumNZ_];
  ColPointers_ = new int[NumMyCols+1];

  // Set pointers into the RowIndices and Values arrays, define ColPointers
  NumIndices = 0;
  for(int i=0; i<NumMyCols; i++) {
    ColPointers_[i] = NumIndices;
    TransIndices_[i] = RowIndices_ + NumIndices;
    TransValues_[i] = RowValues_ + NumIndices;
    NumIndices += TransNumNZ_[i];
  }
  ColPointers_[NumMyCols] = NumIndices;

  Copy();

  int m = A_->NumGlobalRows();
  int n = A_->NumGlobalCols();

  /* Create matrix A in the format expected by SuperLU. */
  dCreate_CompCol_Matrix( &(data_->A), m, n, TotNumNZ_, RowValues_,
	RowIndices_, ColPointers_, SLU_NC, SLU_D, SLU_GE );
  
  /* Create right-hand side matrix B. */
  double * rhs_x;
  double * rhs_b;
  int LDA_x, LDA_b;
  X_->ExtractView( &rhs_x, &LDA_x );
  dCreate_Dense_Matrix( &(data_->X), m, 1, rhs_x, m, SLU_DN, SLU_D, SLU_GE);
  B_->ExtractView( &rhs_b, &LDA_b );
  dCreate_Dense_Matrix( &(data_->B), m, 1, rhs_b, m, SLU_DN, SLU_D, SLU_GE);
  
  perm_r_ = new int[m];
  perm_c_ = new int[n];
  
  etree_ = new int[n];

  ferr_ = new double[1];
  berr_ = new double[1];

  R_ = new double[m];
  C_ = new double[n];
}

Epetra_SLU::~Epetra_SLU()
{
  SLU::Destroy_SuperMatrix_Store( &(data_->A) );
  SLU::Destroy_SuperMatrix_Store( &(data_->B) );
  SLU::Destroy_SuperMatrix_Store( &(data_->X) );

  delete [] TransIndices_;
  delete [] TransValues_;

  delete [] TransNumNZ_;

  delete [] RowIndices_;
  delete [] RowValues_;
  delete [] ColPointers_;

  delete [] perm_r_;
  delete [] perm_c_;

  delete [] etree_;

  delete [] ferr_;
  delete [] berr_;

  delete [] R_;
  delete [] C_;

//  Destroy_CompCol_Matrix(data_->A);
//  SLU::Destroy_SuperNode_Matrix( &(data_->L) );
//  SLU::Destroy_CompCol_Matrix( &(data_->U) );

  delete data_;

  delete AG_;
}

void Epetra_SLU::Copy()
{
  int NumMyCols = AG_->NumMyCols();
  int NumMyEqs = A_->NumMyRows();
  int GlobalMaxNumIndices = AG_->GlobalMaxNumIndices();
  int NumIndices;
  int * xIndices;
  double * xValues;

  for ( int i = 0; i < NumMyCols; i++ )
    TransNumNZ_[i] = 0;

  for ( int i = 0; i < NumMyEqs; i++ )
  {
    assert(A_->ExtractMyRowView( i, NumIndices, xValues, xIndices) == 0 );
    int ii = A_->GRID(i);
    for ( int j = 0; j < NumIndices; j++ )
    {
      int TransRow = xIndices[j];
      int loc = TransNumNZ_[TransRow];
      TransIndices_[TransRow][loc] = ii;
      TransValues_[TransRow][loc] = xValues[j];
      ++TransNumNZ_[TransRow]; // increment counter into current transpose row
    }
  }
}

int Epetra_SLU::Solve( bool Verbose,
                       bool Equil,
                       bool Factor,
                       int perm_type,
                       double pivot_thresh,
                       bool Refact,
                       bool Trans )
{
  Copy();

  using namespace SLU;

  int m = A_->NumGlobalRows();

  int permt = perm_type;
  if( m < 3 ) permt = 0;

  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: use the natural ordering
   *   permc_spec = 1: use minimum degree ordering on structure of A'*A
   *   permc_spec = 2: use minimum degree ordering on structure of A'+A
   */
  if( !count_ ) get_perm_c( permt, &(data_->A), perm_c_ );

  if( Verbose ) cout << "MATRIX COPIED!" << endl;

  int info = 0;

  char fact, trans, refact, equed;

  if( Trans ) trans = 'T';
  else        trans = 'N';

  if( Equil ) fact = 'E';
  else        fact = 'N';

  if( !count_ || !Refact )
  {
    refact = 'N';
    count_ = 0;
  }
  else
  {
    refact = 'Y';
    if( !Factor ) fact = 'F';
  }

  if( Equil ) equed = 'B';
  else        equed = 'N';

  data_->iparam.diag_pivot_thresh = pivot_thresh;

  double rpg, rcond;
 
  if( Verbose ) cout << "TRANS:  " << trans << endl;
  if( Verbose ) cout << "REFACT: " << refact << endl;

  dgssvx( &fact, &trans, &refact, &(data_->A), &(data_->iparam), perm_c_,
	perm_r_, etree_, &equed, R_, C_, &(data_->L), &(data_->U),
	NULL, 0, &(data_->B), &(data_->X), &rpg, &rcond,
	ferr_, berr_, &(data_->mem_usage), &info );

  if( info ) cout << "WARNING: SuperLU returned with error code = " << info << endl;

  if( Verbose )
  {
    cout << "SYSTEM DIRECT SOLVED!" << endl;

    cout << "SuperLU INFO: " << info << "\n\n";
    cout << "  ncols = " << A_->NumGlobalCols() << endl;

    cout << "SuperLU Memory Usage\n";
    cout << "--------------------\n";
    cout << "LU: " << data_->mem_usage.for_lu << endl;
    cout << "Total: " << data_->mem_usage.total_needed << endl;
    cout << "Expands: " << data_->mem_usage.expansions << endl;
    cout << "--------------------\n\n";
  
    if (m<200) dPrint_CompCol_Matrix("A", &(data_->A));
//    if (m<200) dPrint_CompCol_Matrix("U", &(data_->U));
//    if (m<200) dPrint_SuperNode_Matrix("L", &(data_->L));
//  if (m<200) PrintInt10("\nperm_r", m, perm_r);
    if (m<200) dPrint_Dense_Matrix("B", &(data_->B));
    if (m<200) dPrint_Dense_Matrix("X", &(data_->X));
  }

  count_++;

  return 0;
}

