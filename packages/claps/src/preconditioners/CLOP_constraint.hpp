//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef CLOP_CONSTRAINT_HPP
#define CLOP_CONSTRAINT_HPP
#include <stdio.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <math.h>
#include "myzero.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

class CLOP_constraint {
 public:
  CLOP_constraint(const Epetra_CrsMatrix* A_,
		  std::ofstream *fout_,
		  const int print_flag_);
  ~CLOP_constraint();
  int factor();
  void Tran(Epetra_CrsMatrix* & Tran, Epetra_Map* & RowMapMyCon,
	    int* & mycdof, int & nx2, int* & x2_dof, int & nsub_gdofs, 
	    int sub_gdofs[], Epetra_CrsMatrix* & CtT);
 private:
  void determine_if_all_simple(int & simple_flag);
  void determine_nnz_con(int nnz_con[]);
  void get_best_row(int col, int icol, int & iproc);
  int update1(int col, int icol);
  void update2(int col, int icol, int n);
  void update3(int col, int icol, int n, int con_num);
  void update4(int con_num, double a[]);
  int get_ncol(int col, int dir);
  void resize_ncon_me(int n);
  void resize_c_and_s(int n);
  void get_cols_and_sfacs(int col, int icol, int dir);
  void remove_zero(int col2);
  void get_old_info(int col, int & iproc, int & con_num);
  void scale_column(int icol, double sf);
  int find_column(int col);
  void add_columns(int icol1, int col2, double sf);
  void adjust_col_info(int icol1, int icol2);
  void add_new_column(int icol1, double sf, int col2);
  void expand(int* int_temp, int* & a, int n);
  void expand(bool* bool_temp, bool* & a, int n);
  void expand(int** pint_temp, int** & a, int n);
  void expand(double** pdouble_temp, double** & a, int n);

  const Epetra_CrsMatrix* A;
  std::ofstream *fout;
  const int print_flag;
  const Epetra_Comm & Comm;
  Epetra_CrsMatrix *E2T;
  int nrow, ncol, ncol_global, maxcol, *nnzcol, *dimcol, **colrows, blocksize;
  int *ivec, *mycols, *imap, *row_degree, my_best_row, MyPID, max_nnz_con;
  int *cola, *con_row, *con_col, dimsfac, ncon_me, dimcon_me;
  int *NumEntries1, *NumEntries2, **Indices1, **Indices2;
  double **colvals, *dvec, *sfaca, **Values1, **Values2, *max_abs_con;
  bool *row_flag, *con_flag;
};
#endif // CLOP_CONSTRAINT_HPP
  
