#include "CLOP_constraint.hpp"
#include "CRD_utils.hpp"
#include <assert.h>

CLOP_constraint::CLOP_constraint(const Epetra_CrsMatrix* A_) 
  : A(A_), Comm(A->Comm())
{
  int i, j, col, NumEntries, *Indices;
  double *Values;
  //
  // initialize and load data into variables
  //
  blocksize = 5000;
  max_nnz_con = 50;
  nrow = A->NumMyRows();
  ncol = A->NumMyCols();
  ncol_global = A->NumGlobalCols();
  MyPID = Comm.MyPID();
  maxcol = 2*ncol; if (maxcol > ncol_global) maxcol = ncol_global;
  colrows = new int*[maxcol];
  colvals = new double*[maxcol];
  mycols = new int[maxcol];
  nnzcol = new int[maxcol]; myzero(nnzcol, maxcol);
  dimcol = new int[maxcol];
  ivec = new int[nrow];
  dvec = new double[nrow];
  imap = new int[nrow]; memset(imap, -1, nrow*sizeof(int));
  row_degree = new int[nrow];
  row_flag = new bool[nrow]; for (i=0; i<nrow; i++) row_flag[i] = false;
  dimsfac = 20; dimcon_me = 20;
  cola = new int[dimsfac];
  sfaca = new double[dimsfac];
  con_row = new int[dimcon_me];
  con_col = new int[dimcon_me];
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    row_degree[i] = NumEntries;
    for (j=0; j<NumEntries; j++) nnzcol[Indices[j]]++;
  }
  for (i=0; i<ncol; i++) {
    mycols[i] = A->GCID(i);
    dimcol[i] = 2*nnzcol[i]; if (dimcol[i] > nrow) dimcol[i] = nrow;
    colrows[i] = new int[dimcol[i]];
    colvals[i] = new double[dimcol[i]];
  }
  myzero(nnzcol, ncol);
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      col = Indices[j];
      colrows[col][nnzcol[col]] = i;
      colvals[col][nnzcol[col]] = Values[j];
      nnzcol[col]++;
    }
  }
}

CLOP_constraint::~CLOP_constraint()
{
  delete E2T;
}

int CLOP_constraint::factor()
{
  int col, i, iproc, n, icol, con_num, simple_flag;
  //
  // determine if all constraints are simple
  //
  determine_if_all_simple(simple_flag);
  if (MyPID == 0) {
    if (simple_flag != 0) cout << "constraints of type simple " << endl;
    else                  cout << "constraints of type complex" << endl;
  }
  //
  // forward factorization of constraints
  //
  if (simple_flag == 0) {
    ncon_me = 0;
    for (col=0; col<ncol_global; col++) {
      icol = find_column(col);
      get_best_row(col, icol, iproc);
      if (MyPID == iproc) n = update1(col, icol);
      Comm.Broadcast(&n, 1, iproc);
      if (n > dimsfac) resize_c_and_s(n);
      if (MyPID == iproc) update2(col, icol, n);
      Comm.Broadcast(cola, n, iproc);
      Comm.Broadcast(sfaca, n, iproc);
      for (i=0; i<n; i++) if (icol != -1) add_columns(icol, cola[i], sfaca[i]);
      if (MyPID == iproc) for (i=0; i<n; i++) remove_zero(cola[i]);
    }
    //
    // backward factorization of constaints
    //
    for (col=ncol_global-1; col>=0; col--) {
      icol = find_column(col);
      get_old_info(col, iproc, con_num);
      if (MyPID == iproc) n = get_ncol(col, 1);
      Comm.Broadcast(&n, 1, iproc);
      if (n > dimsfac) resize_c_and_s(n);
      if (MyPID == iproc) update3(col, icol, n, con_num);
      Comm.Broadcast(cola, n, iproc);
      Comm.Broadcast(sfaca, n, iproc);
      for (i=0; i<n; i++) if (icol != -1) add_columns(icol, cola[i], sfaca[i]);
      if (MyPID == iproc) for (i=0; i<n; i++) remove_zero(cola[i]);
    }
    //
    // final scaling
    //
    for (col = 0; col<ncol_global; col++) {
      icol = find_column(col);
      get_old_info(col, iproc, con_num);
      if (MyPID == iproc) update4(con_num, sfaca);
      Comm.Broadcast(sfaca, 1, iproc);
      if (icol != -1) scale_column(icol, sfaca[0]);
    }
  }
  //
  // free memory
  //
  delete [] dimcol; delete [] dvec; delete [] imap;
  delete [] row_degree; delete [] cola; delete [] sfaca;
  if (MyPID == 0) cout << "CLOP_constraint::factor completed" << endl;
  return 0;
}

void CLOP_constraint::determine_if_all_simple(int & simple_flag)
{
  int i, j, NumEntries, *Indices, nmpi, nsend, nremain, jpos, *row_con;
  double *Values, *svals, *svals_global, *svals_min, *svals_max;
  svals = new double[ncol]; myzero(svals, ncol);
  svals_global = new double[blocksize]; myzero(svals_global, blocksize);
  svals_min = new double[blocksize];
  svals_max = new double[blocksize];
  row_con = new int[ncol];
  ncon_me = 0;
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    if ((NumEntries == 1) && (fabs(Values[0]) > fabs(svals[Indices[0]]))) {
      svals[Indices[0]] = Values[0];
      row_con[Indices[0]] = i;
    }
  }
  simple_flag = 1;
  nmpi = ncol_global/blocksize;
  nremain = ncol_global - nmpi*blocksize;
  if (nremain > 0) nmpi++;
  for (i=0; i<nmpi; i++) {
    int ibeg = i*blocksize;
    nsend = blocksize;
    if ((i == (nmpi-1)) && (nremain > 0)) nsend = nremain;
    for (j=0; j<ncol; j++) {
      jpos = mycols[j] - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) svals_global[jpos] = svals[j];
    }
    Comm.MinAll(svals_global, svals_min, nsend);
    Comm.MaxAll(svals_global, svals_max, nsend);
    for (j=0; j<nsend; j++) {
      if ((svals_min[j] == 0) && (svals_max[j] == 0)) {
	simple_flag = 0; delete [] svals; delete [] svals_global;
	delete [] svals_min; delete [] svals_max; delete [] row_con; return;
      }
    }
    for (j=0; j<ncol; j++) {
      jpos = mycols[j] - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) {
	svals_global[jpos] = 0;
	double aa = svals_max[jpos];
	if (-svals_min[jpos] > aa) aa = svals_min[jpos];
	if (svals[j] == aa) svals_global[jpos] = MyPID+1;
	svals[j] = 1/aa;
      }
    }
    Comm.MaxAll(svals_global, svals_max, nsend);
    for (j=0; j<ncol; j++) {
      jpos = mycols[j] - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) {
	svals_global[jpos] = 0;
	if (MyPID == (svals_max[jpos]-1)) {
	  if (ncon_me+1 > dimcon_me) resize_ncon_me(ncon_me+1);
	  con_row[ncon_me] = row_con[j];
	  con_col[ncon_me] = j;
	  ncon_me++;
	}
      }
    }
  }
  int simple_flag_min, simple_flag_max;
  Comm.MinAll(&simple_flag, &simple_flag_min, 1);
  Comm.MaxAll(&simple_flag, &simple_flag_max, 1);
  assert(simple_flag == simple_flag_min);
  assert(simple_flag == simple_flag_max);
  for (j=0; j<ncol; j++) {
    for (i=0; i<nnzcol[j]; i++) colvals[j][i] *= svals[j];
  }
  delete [] svals; delete [] svals_global; delete [] svals_min; delete [] svals_max;
  delete [] row_con;
}

int CLOP_constraint::update1(int col, int icol)
{
  int n;
  assert (icol != -1);
  n = get_ncol(col, 0);
  if (ncon_me+1 > dimcon_me) resize_ncon_me(ncon_me+1);
  con_row[ncon_me] = my_best_row;
  con_col[ncon_me] = icol;
  return n;
}

void CLOP_constraint::update2(int col, int icol, int n)
{
  get_cols_and_sfacs(col, icol, 0);
  ncon_me++;
}

void CLOP_constraint::update3(int col, int icol, int n, int con_num)
{
  get_cols_and_sfacs(col, icol, 1);
}

void CLOP_constraint::update4(int con_num, double a[])
{
  int i, row, icol, row2;
  row = con_row[con_num];
  icol = con_col[con_num];
  for (i=0; i<nnzcol[icol]; i++) {
    row2 = colrows[icol][i];
    if (row2 == row) {
      a[0] = 1.0/colvals[icol][i];
      break;
    }
  }
}

void CLOP_constraint::Tran(Epetra_CrsMatrix* & Tran, Epetra_Map* & RowMapMyCon,
			   int* & mycdof, int & nx2, int* & x2_dof, int & nsub_gdofs, 
			   int sub_gdofs[], Epetra_CrsMatrix* & CtT)
{
  int i, j, NumEntries, row, ncolT, gdof, nc1, nc2, nc2_all, *Indices, ierr;
  double *Values;
  //
  // flag constrained dofs and determine Epetra_Map for "owned" constraints
  //
  for (i=0; i<nrow; i++) row_flag[i] = false;
  int *mycons = new int[ncon_me];
  mycdof = new int[ncon_me];
  for (i=0; i<ncon_me; i++) {
    mycons[i] = mycols[con_col[i]];
    mycdof[i] = con_row[i];
    row_flag[con_row[i]] = true;
  }
  RowMapMyCon = new Epetra_Map(ncol_global, ncon_me, mycons, 0, Comm);
  //
  // form transpose of transformed constraint matrix
  //
  myzero(ivec, nrow);
  for (i=0; i<ncon_me; i++) ivec[con_row[i]] = 1;
  for (i=0; i<ncol; i++) {
    for (j=0; j<nnzcol[i]; j++) {
      row = colrows[i][j];
      if (row_flag[row] == false) ivec[row]++;
    }
  }
  int *val2 = new int[nrow+1]; val2[0] = 0;
  for (i=0; i<nrow; i++) val2[i+1] = val2[i] + ivec[i];
  int *cval1 = new int[val2[nrow]];
  double *val1 = new double[val2[nrow]];
  myzero(ivec, nrow);
  for (i=0; i<ncon_me; i++) {
    val1[ val2[con_row[i]]] = 1;
    cval1[val2[con_row[i]]] = con_col[i];
    ivec[con_row[i]] = 1;
  }
  double con_tol(1e-10);
  for (i=0; i<ncol; i++) {
    for (j=0; j<nnzcol[i]; j++) {
      row = colrows[i][j];
      if (row_flag[row] == false) {
	if (fabs(colvals[i][j]) > con_tol) {
	  val1[ val2[row] + ivec[row]] = colvals[i][j];
	  cval1[val2[row] + ivec[row]] = i;
	  ivec[row]++;
	}
      }
    }
  }
  Epetra_Map NewColMap(-1, ncol, mycols, 0, Comm);
  CtT = new Epetra_CrsMatrix(Copy, A->RowMap(), NewColMap, ivec);
  for (i=0; i<nrow; i++) {
    ierr = CtT->InsertMyValues(i, ivec[i], &val1[val2[i]], &cval1[val2[i]]);
    assert (ierr == 0);
  }
  CtT->FillComplete(A->DomainMap(), A->RangeMap());
  /*
  char filename[101];
  sprintf(filename,"%s%d","CtT", MyPID);
  CRD_utils::Epetra_datfile(CtT, filename);
  sprintf(filename,"%s%d","nnzcol", MyPID);
  CRD_utils::Epetra_datfile(nnzcol, ncol, filename);
  */
  delete [] ivec; delete [] val1; delete [] val2; delete [] cval1;
  //
  // construct Epetra_CrsMatrix for on-processor contributions to constraint matrix
  //
  Epetra_CrsMatrix ConMat_Loc(View, NewColMap, A->RowMap(), 0);
  for (i=0; i<ncol; i++) {
    NumEntries = 0;
    for (j=0; j<nnzcol[i]; j++) {
      row = colrows[i][j];
      if (row_flag[row] == false) {
	if (fabs(colvals[i][j]) > con_tol) {
	  colrows[i][NumEntries] = colrows[i][j];
	  colvals[i][NumEntries] = colvals[i][j];
	  NumEntries++;
	}
      }
    }
    ierr = ConMat_Loc.InsertMyValues(i, NumEntries, colvals[i], colrows[i]);
    assert (ierr == 0);
  }
  ConMat_Loc.FillComplete(A->RowMap(), *RowMapMyCon);
  //
  // determine number of nonzeros (global) for each constraint
  //  nnz_con[i]  = number of nonzeros in global constraint mycols[i]
  //
  int *nnz_con = new int[ncol];
  determine_nnz_con(nnz_con);
  //
  // determine type 1 constraints and adjust row_flag to only flag 
  // type 1 constrained dofs
  //
  for (i=0; i<nrow; i++) row_flag[i] = false;
  nc1 = 0;
  for (i=0; i<ncon_me; i++) {
    if (nnz_con[con_col[i]] <= max_nnz_con) {
      mycons[nc1] = mycols[con_col[i]];
      row_flag[con_row[i]] = true;
      nc1++;
    }
  }
  //
  // determine non-constrained (type 1) dofs in sub_gdofs array
  //
  int *row_active = new int[nrow]; myzero(row_active, nrow);
  int *row_active_sub = new int[nsub_gdofs]; myzero(row_active_sub, nsub_gdofs);
  for (i=0; i<nrow; i++) if (row_flag[i] == true) row_active[i] = 1;
  Epetra_Map RowMap_sub(-1, nsub_gdofs, sub_gdofs, 0, Comm);
  Epetra_IntVector RA(View, A->RowMap(), row_active);
  Epetra_IntVector RA_sub(View, RowMap_sub, row_active_sub);
  Epetra_Import Importer(RowMap_sub, A->RowMap());
  RA_sub.Import(RA, Importer, Insert);
  int nsub_gdofs_orig = nsub_gdofs; nsub_gdofs = 0;
  for (i=0; i<nsub_gdofs_orig; i++) {
    if (row_active_sub[i] == 0) {
      sub_gdofs[nsub_gdofs] = sub_gdofs[i];
      nsub_gdofs++;
    }
  }
  delete [] row_active;
  delete [] row_active_sub;
  //
  // gather contributions to type 1 constraints
  //
  Epetra_Map RowMapMyCon1(-1, nc1, mycons, 0, Comm);
  Epetra_CrsMatrix ConMat1(Copy, RowMapMyCon1, 0);
  Epetra_Export Exporter(ConMat_Loc.RowMap(), RowMapMyCon1);
  ConMat1.Export(ConMat_Loc, Exporter, Insert);
  ConMat1.FillComplete(A->RowMap(), RowMapMyCon1);
  //  cout << ConMat1 << endl;
  //
  // form transformation matrix for type 1 constraints
  //
  int ncol1 = nrow - nc1;
  int ncol2 = ConMat1.NumMyCols();
  int ncol_all = ncol1 + ncol2;
  int *gcol_all = new int[ncol_all];
  ncol_all = 0;
  for (i=0; i<nrow; i++) {
    if (row_flag[i] == false) {
      gcol_all[ncol_all] = A->GRID(i);
      ncol_all++;
    }
  }
  int ndof_u = ncol_all;
  for (i=0; i<ncol2; i++) {
    gcol_all[ncol_all] = ConMat1.GCID(i);
    ncol_all++;
  }
  int *dof_u = new int[ndof_u]; ndof_u = 0;
  for (i=0; i<nrow; i++) {
    if (row_flag[i] == false) {
      dof_u[ndof_u] = A->GRID(i);
      ndof_u++;
    }
  }
  CRD_utils::sort_and_cull(gcol_all, ncol_all, ncolT);
  Epetra_Map ColMapT(-1, ncolT, gcol_all, 0, Comm);
  Tran = new Epetra_CrsMatrix(Copy, A->RowMap(), ColMapT, 0);
  double one(1.0);
  for (i=0; i<nrow; i++) {
    if (row_flag[i] == false) {
      gdof = A->GRID(i);
      int loc_col = CRD_utils::find_index(gcol_all, ncolT, gdof);
      assert (loc_col != -1);
      ierr = Tran->InsertMyValues(i, 1, &one, &loc_col);
      assert (ierr == 0);
    }
  }
  nc1 = 0;
  int *Iarray = new int[max_nnz_con];
  double *Darray = new double[max_nnz_con];
  for (i=0; i<ncon_me; i++) {
    if (nnz_con[con_col[i]] <= max_nnz_con) {
      ConMat1.ExtractMyRowView(nc1, NumEntries, Values, Indices);
      assert (NumEntries <= max_nnz_con);
      for (j=0; j<NumEntries; j++) {
	gdof = ConMat1.GCID(Indices[j]);
	int loc_col = CRD_utils::find_index(gcol_all, ncolT, gdof);
	assert (loc_col != -1);
	Iarray[j] = loc_col;
	Darray[j] = -Values[j];
      }
      ierr = Tran->InsertMyValues(con_row[i], NumEntries, Darray, Iarray);
      assert (ierr == 0);
      nc1++;
    }
  }
  delete [] Iarray; delete [] Darray;
  Epetra_Map RowMapu(-1, ndof_u, dof_u, 0, Comm);
  delete [] dof_u; delete [] gcol_all; delete [] row_flag; 
  Tran->FillComplete(RowMapu, A->RowMap());
  /*
  sprintf(filename,"%s%d","Tran", MyPID);
  CRD_utils::Epetra_datfile(Tran, filename);
  */
  //  cout << *Tran << endl;
  //
  // construct Epetra_CrsMatrix for on-processor contribution to transpose of
  // type 2 constraints (E_2^T in notes)
  //
  nc2 = 0;
  for (i=0; i<ncon_me; i++) {
    if (nnz_con[con_col[i]] > max_nnz_con) {
      mycons[nc2] = mycols[con_col[i]];
      nc2++;
    }
  }
  x2_dof = new int[nc2]; nc2 = 0;
  for (i=0; i<ncon_me; i++) {
    if (nnz_con[con_col[i]] > max_nnz_con) {
      x2_dof[nc2] = con_row[i];
      nc2++;
    }
  }
  nx2 = nc2;
  nc2_all = 0;
  for (i=0; i<ncol; i++) if (nnz_con[i] > max_nnz_con) nc2_all++;
  int *col_t2 = new int[nc2_all]; nc2_all = 0;
  for (i=0; i<ncol; i++) {
    if (nnz_con[i] > max_nnz_con) {
      col_t2[nc2_all] = mycols[i];
      nc2_all++;
    }
  }
  Epetra_Map RowMapMyCon2(-1, nc2, mycons, 0, Comm);
  Epetra_Map RowMapCon2(-1, nc2_all, col_t2, 0, Comm);
  E2T = new Epetra_CrsMatrix(Copy, A->RowMap(), RowMapCon2, 0);
  delete [] col_t2; delete [] mycons;
  nc2_all = 0;
  for (i=0; i<ncol; i++) {
    if (nnz_con[i] > max_nnz_con) {
      ConMat_Loc.ExtractMyRowView(i, NumEntries, Values, Indices);
      for (j=0; j<NumEntries; j++) {
	ierr = E2T->InsertMyValues(Indices[j], 1, &Values[j], &nc2_all);
	assert (ierr == 0);
      }
      nc2_all++;
    }
  }
  for (i=0; i<ncol; i++) {
    delete [] colrows[i];
    delete [] colvals[i];
  }
  delete [] colrows; delete [] colvals; delete [] mycols; delete [] nnzcol;
  delete [] con_row; delete [] con_col; delete [] nnz_con;
}

void CLOP_constraint::determine_nnz_con(int nnz_con[])
{
  int i, j, nmpi, nremain, nsend, ibeg, jpos;
  int *nnz_sum = new int[blocksize];
  int *nnz_loc = new int[blocksize]; myzero(nnz_loc, blocksize);
  nmpi = ncol_global/blocksize;
  nremain = ncol_global - nmpi*blocksize;
  if (nremain > 0) nmpi++;
  for (i=0; i<nmpi; i++) {
    ibeg = i*blocksize;
    nsend = blocksize;
    if ((i == (nmpi-1)) && (nremain > 0)) nsend = nremain;
    for (j=0; j<ncol; j++) {
      jpos = mycols[j] - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) nnz_loc[jpos] = nnzcol[j];
    }
    Comm.SumAll(nnz_loc, nnz_sum, nsend);
    for (j=0; j<ncol; j++) {
      jpos = mycols[j] - ibeg;
      if ((jpos >= 0) && (jpos < nsend)) {
	nnz_loc[jpos] = 0;
	nnz_con[j] = nnz_sum[jpos];
      }
    }
  }
  delete [] nnz_sum; delete [] nnz_loc;
}

void CLOP_constraint::get_best_row(int col, int icol, int & iproc)
{
  int i, row;
  double max_value(0), colv, max_value_g, tol(1e-3);
  //
  // find largest entry in column
  //
  if (icol != -1) {
    for (i=0; i<nnzcol[icol]; i++) {
      colv = colvals[icol][i];
      if (fabs(colv) > max_value) max_value = fabs(colv);
    }
  }
  Comm.MaxAll(&max_value, &max_value_g, 1);
  assert (max_value_g > 0);
  //
  // determine lowest possible degree
  //
  int max_degree(7777777); int best_degree = max_degree;
  int best_degree_g, iflag(0);
  if (icol != -1) {
    for (i=0; i<nnzcol[icol]; i++) {
      colv = colvals[icol][i];
      if (fabs(colv) > tol*max_value_g) {
	row = colrows[icol][i];
	if ((row_degree[row] < best_degree) &&
	    (row_flag[row] == false)) best_degree = row_degree[row];
      }
    }
  }
  Comm.MinAll(&best_degree, &best_degree_g, 1);
  if (best_degree_g == max_degree) {
    iflag = 1;
    if (icol != -1) {
      for (i=0; i<nnzcol[icol]; i++) {
	colv = colvals[icol][i];
	if (fabs(colv) > tol*max_value_g) {
	  row = colrows[icol][i];
	  if (row_degree[row] < best_degree) best_degree = row_degree[row];
	}
      }
    }
    Comm.MinAll(&best_degree, &best_degree_g, 1);
  }
  assert (best_degree_g < max_degree);
  //
  // determine pivot row
  //
  max_value = 0; my_best_row = -1;
  if (icol != -1) {
    for (i=0; i<nnzcol[icol]; i++) {
      colv = colvals[icol][i];
      if (fabs(colv) > tol*max_value_g) {
	row  = colrows[icol][i];
	if ((row_degree[row] == best_degree) && (fabs(colv) > max_value)) {
	  if ((iflag == 1) || (row_flag[row] == false)) {
	    max_value = fabs(colv);
	    my_best_row = row;
	  }
	}
      }
    }
  }
  Comm.MaxAll(&max_value, &max_value_g, 1);
  if (max_value == max_value_g) max_value = (double) MyPID;
  else max_value = -1;
  Comm.MaxAll(&max_value, &max_value_g, 1);
  iproc = -1;
  if (max_value == max_value_g) iproc = MyPID;
  int maxproc;
  Comm.MaxAll(&iproc, &maxproc, 1);
  iproc = maxproc;
  assert (iproc >= 0);
}

int CLOP_constraint::get_ncol(int col, int dir)
{
  //
  // for column col determine number of other columns (to right or left)
  // with a nonzero entry in row my_best_row
  //
  int i, j, n(0), col2;
  for (i=0; i<ncol; i++) {
    col2 = mycols[i];
    if (((col2 > col) && (dir == 0)) || ((col2 < col) && (dir == 1))) {
      for (j=0; j<nnzcol[i]; j++) {
	if (colrows[i][j] == my_best_row) n++;
      }
    }
  }
  return n;
}
   
void CLOP_constraint::resize_ncon_me(int n)
{
  int i;
  dimcon_me = 2*n; 
  if (dimcon_me > ncol_global) dimcon_me = ncol_global;
  int *con_row_temp = new int[ncon_me];
  int *con_col_temp = new int[ncon_me];
  for (i=0; i<ncon_me; i++) {
    con_row_temp[i] = con_row[i]; con_col_temp[i] = con_col[i];
  }
  delete [] con_row; delete [] con_col; 
  con_row = new int[dimcon_me]; con_col = new int[dimcon_me];
  for (i=0; i<ncon_me; i++) {
    con_row[i] = con_row_temp[i]; con_col[i] = con_col_temp[i];
  }
  delete [] con_row_temp; delete [] con_col_temp;
}

void CLOP_constraint::resize_c_and_s(int n)
{
  dimsfac = 2*n; 
  if (dimsfac > ncol_global) dimsfac = ncol_global;
  delete [] cola; delete [] sfaca;
  cola = new int[dimsfac]; sfaca = new double[dimsfac];
}

void CLOP_constraint::get_cols_and_sfacs(int col, int icol, int dir)
{
  int i, j, n(0), row, col2;
  double pivot_val(0.0);
  for (j=0; j<nnzcol[icol]; j++) {
    row = colrows[icol][j];
    if (row == my_best_row) pivot_val = colvals[icol][j];
  }
  assert (pivot_val != 0);
  for (i=0; i<ncol; i++) {
    col2 = mycols[i];
    if (((col2 > col) && (dir == 0)) || ((col2 < col) && (dir == 1))) {
      for (j=0; j<nnzcol[i]; j++) {
	if (colrows[i][j] == my_best_row) {
	  cola[n] = col2;
	  sfaca[n] = -colvals[i][j]/pivot_val;
	  n++;
	}
      }
    }
  }
}

void CLOP_constraint::remove_zero(int col2)
{
  int i, j, nnz(0), row;
  i = find_column(col2);
  for (j=0; j<nnzcol[i]; j++) {
    if (colrows[i][j] != my_best_row) {
      colrows[i][nnz] = colrows[i][j];
      colvals[i][nnz] = colvals[i][j];
      nnz++;
    }
  }
  nnzcol[i] = nnz;
}   

void CLOP_constraint::get_old_info(int col, int & iproc, int & con_num)
{
  int i, maxproc;
  iproc = -1;
  for (i=0; i<ncon_me; i++) {
    if (mycols[con_col[i]] == col) { 
      iproc = MyPID;
      my_best_row = con_row[i];
      con_num = i;
    }
  }
  Comm.MaxAll(&iproc, &maxproc, 1);
  iproc = maxproc;
  assert (iproc >= 0);
}

void CLOP_constraint::scale_column(int icol, double sf)
{
  int i;
  for (i=0; i<nnzcol[icol]; i++) colvals[icol][i] *= sf;
}

int CLOP_constraint::find_column(int col)
{
  //
  // inefficient search to find column col in mycols array
  //
  for (int i=0; i<ncol; i++) if (mycols[i] == col) return i;
  return -1;
}

void CLOP_constraint::add_columns(int icol1, int col2, double sf)
{
  //
  // col2 values = sf * col1 values + col2 values
  //
  int i, row1, irow2, icol2;
  icol2 = find_column(col2);
  if (icol2 == -1) add_new_column(icol1, sf, col2);
  else {
    for (i=0; i<nnzcol[icol2]; i++) imap[colrows[icol2][i]] = i;
    adjust_col_info(icol1, icol2);
    for (i=0; i<nnzcol[icol1]; i++) {
      row1 = colrows[icol1][i];
      irow2 = imap[row1];
      assert(irow2 != -1);
      colvals[icol2][irow2] += sf*colvals[icol1][i];
    }
    for (i=0; i<nnzcol[icol2]; i++) imap[colrows[icol2][i]] = -1;
  }
}

void CLOP_constraint::adjust_col_info(int icol1, int icol2)
{
  int i, row, nadd(0);
  for (i=0; i<nnzcol[icol1]; i++) {
    row = colrows[icol1][i];
    if (imap[row] == -1) nadd++;
  }
  int nnew = nnzcol[icol2] + nadd;
  if (nnew > dimcol[icol2]) {
    dimcol[icol2] = 2*nnew; if (dimcol[icol2] > nrow) dimcol[icol2] = nrow;
    memcpy(dvec, colvals[icol2], nnzcol[icol2]*sizeof(double));
    memcpy(ivec, colrows[icol2], nnzcol[icol2]*sizeof(int));
    delete [] colvals[icol2];
    delete [] colrows[icol2];
    colvals[icol2] = new double[dimcol[icol2]];
    colrows[icol2] = new int[dimcol[icol2]];
    memcpy(colvals[icol2], dvec, nnzcol[icol2]*sizeof(double));
    memcpy(colrows[icol2], ivec, nnzcol[icol2]*sizeof(int));
  }
  for (i=0; i<nnzcol[icol1]; i++) {
    row = colrows[icol1][i];
    if (imap[row] == -1) {
      colrows[icol2][nnzcol[icol2]] = row;
      colvals[icol2][nnzcol[icol2]] = 0;
      imap[row] = nnzcol[icol2];
      nnzcol[icol2]++;
    }
  }
}

void CLOP_constraint::add_new_column(int icol1, double sf, int col2)
{
  int i, nnew;
  nnew = ncol + 1;
  if (nnew > maxcol) {
    maxcol = 2*nnew; if (maxcol > ncol_global) maxcol = ncol_global;
    int *int_temp = new int[ncol];
    int **pint_temp = new int*[ncol];
    double **pdouble_temp = new double*[ncol];
    expand(int_temp, dimcol, maxcol);
    expand(int_temp, nnzcol, maxcol);
    expand(int_temp, mycols, maxcol);
    expand(pint_temp, colrows, maxcol);
    expand(pdouble_temp, colvals, maxcol);
    delete [] int_temp;
    delete [] pint_temp;
    delete [] pdouble_temp;
  }
  mycols[ncol] = col2;
  nnzcol[ncol] = nnzcol[icol1];
  dimcol[ncol] = 2*nnzcol[icol1]; if (dimcol[ncol] > nrow) dimcol[ncol] = nrow;
  colvals[ncol] = new double[dimcol[ncol]];
  colrows[ncol] = new int[dimcol[ncol]];
  for (i=0; i<nnzcol[icol1]; i++) {
    colvals[ncol][i] = sf*colvals[icol1][i];
    colrows[ncol][i] = colrows[icol1][i];
  }
  ncol++;
}

void CLOP_constraint::expand(int* int_temp, int* & a, int n)
{
  int i;
  for (i=0; i<ncol; i++) int_temp[i] = a[i];
  delete [] a; 
  a = new int[n];
  for (i=0; i<ncol; i++) a[i] = int_temp[i];
}

void CLOP_constraint::expand(int** pint_temp, int** & a, int n)
{
  int i;
  for (i=0; i<ncol; i++) pint_temp[i] = a[i];
  delete [] a; 
  a = new int*[n];
  for (i=0; i<ncol; i++) a[i] = pint_temp[i];
}

void CLOP_constraint::expand(double** pdouble_temp, double** & a, int n)
{
  int i;
  for (i=0; i<ncol; i++) pdouble_temp[i] = a[i];
  delete [] a; 
  a = new double*[n];
  for (i=0; i<ncol; i++) a[i] = pdouble_temp[i];
}

