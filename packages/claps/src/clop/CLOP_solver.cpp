#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <assert.h>
#include <climits>
#include "CLOP_solver.hpp"
#include "myzero.hpp"
#include <algorithm>


extern "C"{
  void metis_partgraphrecursive(int *n, int *xadj, int *adjncy, int *vwgt,
				int *adjwgt, int *wgtflag, int *numflag, 
				int *nparts, int *options, int *edgecut, 
				int *part);
  void metis_partgraphkway(int *n, int *xadj, int *adjncy, int *vwgt,
			   int *adjwgt, int *wgtflag, int *numflag, 
			   int *nparts, int *options, int *edgecut, int *part);
  void dstev_(char* N, int* i, double Dtri[], double Etri[], double* Z, 
	      int* one, double* WORK, int* INFO, long lengthN); 
}

CLOP_solver::CLOP_solver(const Epetra_CrsMatrix* AStandard_,
			 const Epetra_IntVector* LDStandard_,
			 const Epetra_MultiVector* CStandard_,
			 const Epetra_Map* SubMap_,
			 const Epetra_CrsMatrix* ConStandard_,
			 const Epetra_IntVector* GNStandard_,
			 const double* clop_params_)
  : AStandard(AStandard_), LDStandard(LDStandard_), 
    CStandard(CStandard_), SubMap(SubMap_), ConStandard(ConStandard_), 
    GNStandard(GNStandard_), clop_params(clop_params_), Comm(AStandard->Comm())
{
  overlap      = int(clop_params[0]); 
  solver_tol   = clop_params[1];
  maxiter      = int(clop_params[2]);
  max_orthog   = int(clop_params[3]);
  atype        = int(clop_params[4]);
  ndim         = int(clop_params[5]);
  local_solver = int(clop_params[6]);
  prt_debug    = int(clop_params[7]);
  prt_summary  = int(clop_params[8]);
  double starttime, endtime;
  //  const Epetra_MpiComm &empicomm = 
  //               dynamic_cast<const Epetra_MpiComm &>(Comm);
  //  mpicomm = empicomm.Comm();
  Comm.Barrier();
  MyPID = Comm.MyPID();
  if (MyPID == 0) starttime = MPI_Wtime();
  NumProc = Comm.NumProc();
  ndof_Standard = AStandard->NumMyRows();
  ndof_global = AStandard->NumGlobalRows();
  print_flag = -1;
  if (MyPID == 0) {
    print_flag = prt_summary + 10*prt_debug;
    fout.open("CLOP_solver.data");
    fout << "----------------- CLOP solver summary information "
	 << "-----------------" << endl;
  }
  if (print_flag > 0) fout << "number of global dofs        = " << ndof_global 
			   << endl;
  zero_pointers();
  //
  // process constraints
  //
  process_constraints();
  //
  // construct AOverlap, LDOverlap, and COverlap
  //
  construct_Overlap();
  //  cout << *AOverlap << endl;
  //  cout << *LDOverlap << endl;
  //  cout << *COverlap << endl;
  construct_subdomains();
  initialize_subdomains();
  //
  // calculate coarse stiffness matrix
  //
  calculate_coarse_stiff();
  //  cout << *Kc << endl;
  //
  // gather coarse stiffness matrix
  //
  gather_coarse_stiff();
  //
  // factor coarse stiffness matrix
  //
  factor_coarse_stiff();
  //
  // initialize solver
  //
  solve_init();
  Comm.Barrier();
  if (print_flag > 0) {
    endtime = MPI_Wtime();
    fout << "elapsed time for clop solver init  = " 
	 << endtime-starttime << " seconds" << endl;
    fout.close();
  }
}

CLOP_solver::~CLOP_solver()
{
  delete AOverlap;
  delete LDOverlap; delete COverlap; delete [] dofpart1; delete [] dofpart2;
  delete [] csdima; delete [] Asub; delete [] cs_local; delete [] xcent;
  delete [] ycent; delete [] zcent; delete [] cdof_proc; delete PhiT;
  delete Phi; delete Kc; delete Kc_gathered; delete Importer_coarse;
  delete ImporterST2O; delete ExporterO2ST; delete Kc_fac; delete GNOverlap;
  delete [] r_overlap; delete [] z_overlap;
  delete [] rhs_coarse; delete [] sol_coarse;
  delete [] temp_coarse; delete PhiTr;
  delete PhiTr_gathered; delete CSol_gathered; delete [] rcurra;
  delete RowMap_coarse; delete rSt_red; delete ApSt_red; 
  delete [] rhs_work; delete [] sol_work; delete [] tmp_work;
  delete rOverlap; delete zOverlap; delete zSt_red; delete pSt_red;
  delete vSt_red; delete [] rhoa; delete [] betaa; delete [] pApa;
  delete [] Etri; delete [] Dtri; delete [] econa; delete ConError;
  delete Tran; delete ASt_red_keep; delete uSt_red; delete vStand; 
  delete wStand; delete [] x2_dof; delete Lambda_local; delete Lambda; 
  delete RowMapMyCon; delete [] lambda_local; delete [] mycdof; 
  delete [] lambda; delete Exporter_lam; delete CtT; delete AP_matrix; 
  delete P_matrix; delete [] pAp_vec; delete [] ortho_vec;
  delete [] ortho_sum; delete [] VV; delete [] HH; delete [] RR;
  delete [] zz; delete [] cc; delete [] ss; delete [] norms;
  delete [] gmres_vec; delete [] gmres_sum; delete gSt_red; delete [] PAP;
  delete [] PAP_sum; delete [] IPIV; delete [] PAP_store;
}

void CLOP_solver::zero_pointers()
{
  AOverlap = 0; LDOverlap = 0; GNOverlap = 0;
  COverlap = 0; dofpart1 = 0; dofpart2 = 0; imap = 0; count1 = 0; csdima = 0;
  Asub = 0; cs_local = 0; xcent = 0; ycent = 0; zcent = 0;
  cdof_proc = 0; PhiT = 0; Phi = 0; Kc = 0; Kc_gathered = 0; 
  Importer_coarse = 0; ImporterST2O = 0; ExporterO2ST = 0; Kc_fac = 0;
  rhs_coarse = 0; sol_coarse = 0; temp_coarse = 0;
  PhiTr = 0; PhiTr_gathered = 0;
  CSol_gathered = 0; rcurra = 0; RowMap_coarse = 0; rSt_red = 0; ApSt_red = 0;
  r_overlap = 0; rhs_work = 0; sol_work = 0; tmp_work = 0; z_overlap = 0;
  rOverlap = 0; zOverlap = 0; zSt_red = 0; pSt_red = 0; vSt_red = 0; rhoa = 0;
  betaa = 0; pApa = 0; Etri = 0; Dtri = 0; econa = 0; ConError = 0;
  Tran = 0; ASt_red_keep = 0; sub_gdofs = 0; uSt_red = 0; vStand = 0; 
  wStand = 0; x2_dof = 0; mycdof = 0; RowMapMyCon = 0; Lambda_local = 0; 
  Lambda = 0; lambda_local = 0; lambda = 0; Exporter_lam = 0;
  CtT = 0; AP_matrix = 0; P_matrix = 0; pAp_vec = 0; ortho_vec = 0;
  ortho_sum = 0; VV = 0; HH = 0; RR = 0; zz = 0; cc = 0; ss = 0; norms = 0;
  gmres_sum = 0; gmres_vec = 0; gSt_red = 0; PAP = 0; PAP_sum = 0; IPIV = 0;
  PAP_store = 0;
}

void CLOP_solver::process_constraints()
{
  int i, j, nrow, ncol, col, NumEntries, *Indices, ierr;
  double *Values;
  //
  // set pointers and return if there are no constraint equations
  //
  if (!SubMap) {
    nsub_gdofs = ndof_Standard;
    sub_gdofs = new int[nsub_gdofs];
    for (i=0; i<ndof_Standard; i++) sub_gdofs[i] = AStandard->GRID(i);
  }
  else {
    nsub_gdofs = SubMap->NumMyPoints();
    sub_gdofs = new int[nsub_gdofs];
    SubMap->MyGlobalElements(sub_gdofs);
  }
  ncon_global = 0;
  if (ConStandard) ncon_global = ConStandard->NumGlobalCols();
  if (print_flag > 0) fout << "number of global constraints = " 
			   << ncon_global << endl;
  if (ncon_global == 0) {
    ASt_red = AStandard;
    return;
  }
  /*
  char filename[101];
  sprintf(filename,"%s%d","CLOP_constraint", MyPID);
  CRD_utils::Epetra_datfile(ConStandard, filename);
  */
  //
  // transform constraints to standard form
  //
  ierr = transform_constraints();
  //
  // calculate reduced stiffness matrix
  //
  Epetra_CrsMatrix *TranT;
  construct_transpose(Tran, TranT);
  Epetra_CrsMatrix *KTran;
  EpetraExtCD::MatrixMatrix::Multiply(*AStandard, *Tran, KTran);
  EpetraExtCD::MatrixMatrix::Multiply(*TranT, *KTran, ASt_red_keep);
  //  cout << *ASt_red_keep << endl;
  ASt_red = ASt_red_keep;

  uSt_red = new Epetra_Vector(ASt_red->RowMap());
  vStand = new Epetra_Vector(AStandard->RowMap());
  uSt_red->Random();
  Tran->Multiply(false, *uSt_red, *vStand);
  ConStandard->Multiply(true, *vStand, *Lambda);
  double inf_norm_error; Lambda->NormInf(&inf_norm_error);
  double con_error_norm = inf_norm_error/ConStandard->NormInf();
  int max_nnz_row = Tran->GlobalMaxNumEntries();
  double inf_norm_Tran = Tran->NormInf();
  if (MyPID == 0) {
    double nnz_before = AStandard->NumGlobalNonzeros();
    double nnz_after  =   ASt_red->NumGlobalNonzeros();
    if (print_flag > 0) {
      fout << "constraint data ------------------------------" << endl;
      fout << "ratio of nnzs after static condensation = "
	   << nnz_after/nnz_before << endl;
      fout << "normalized constraint error check       = " << con_error_norm 
	   << endl;
      fout << "maximum nnz in any row of T matrix      = " << max_nnz_row 
	   << endl;
      fout << "infinity norm of T matrix               = " << inf_norm_Tran 
	   << endl;
    }
  }
  assert(con_error_norm < 1e-10);
  delete KTran;
  delete TranT;
}

int CLOP_solver::transform_constraints()
{
  int flag(0);
  CLOP_constraint *A;
  A = new CLOP_constraint(ConStandard, &fout, print_flag);
  flag = A->factor();
  A->Tran(Tran, RowMapMyCon, mycdof, nx2, x2_dof, nsub_gdofs, sub_gdofs, CtT);
  nmycon = RowMapMyCon->NumMyElements();
  int nn = ConStandard->DomainMap().NumMyPoints();
  lambda_local = new double[nmycon];
  lambda  = new double[nn];
  Lambda_local = new Epetra_Vector(View, *RowMapMyCon, lambda_local);
  Lambda       = new Epetra_Vector(View, ConStandard->DomainMap(), lambda);
  ConError     = new Epetra_Vector(ConStandard->DomainMap());
  Exporter_lam = new Epetra_Export(*RowMapMyCon, ConStandard->DomainMap());
  Comm.SumAll(&nx2, &nx2_global, 1);
  assert(nx2_global == 0);
  delete A;
  return flag;
}

void CLOP_solver::construct_transpose(Epetra_CrsMatrix* & A, 
				      Epetra_CrsMatrix* & AT)
{
  int i, j, k, nrow, ncol, nnz, row, NumEntries, *Indices;
  double *Values;
  nrow = A->RowMap().NumMyElements();
  ncol = A->ColMap().NumMyElements();
  nnz = A->NumMyNonzeros();
  Epetra_CrsMatrix B(View, A->ColMap(), A->RowMap(), 0);
  int *count = new int[ncol]; myzero(count, ncol);
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) count[Indices[j]]++;
  }
  int *rowbeg = new int[ncol+1]; rowbeg[0] = 0;
  int *colidx = new int[nnz];
  double *Bval = new double[nnz];
  for (i=0; i<ncol; i++) rowbeg[i+1] = rowbeg[i] + count[i];
  myzero(count, ncol);
  for (i=0; i<nrow; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      row = Indices[j];
      k = rowbeg[row] + count[row];
      colidx[k] = i;
      Bval[  k] = Values[j];
      count[row]++;
    }
  }
  for (i=0; i<ncol; i++) {
    NumEntries = rowbeg[i+1] - rowbeg[i];
    k = rowbeg[i];
    B.InsertMyValues(i, NumEntries, &Bval[k], &colidx[k]);
  }
  B.FillComplete(A->RangeMap(), A->DomainMap());
  AT = new Epetra_CrsMatrix(Copy, A->DomainMap(), 0);
  Epetra_Export Exporter(B.RowMap(), A->DomainMap());
  AT->Export(B, Exporter, Add);
  AT->FillComplete(A->RangeMap(), A->DomainMap());
  delete [] count; delete [] rowbeg; delete [] colidx; delete [] Bval;
}

void CLOP_solver::construct_Overlap()
{
  if (print_flag > 9) fout << "in construct_Overlap " << endl;
  ndof_global_red = ASt_red->NumGlobalRows();
  int i;
  assert (overlap >= 0);
  Epetra_CrsGraph *OverlapGraph;
  Epetra_Import *Importer;
  Epetra_BlockMap *RowMap;
  RowMap = new Epetra_BlockMap(-1, nsub_gdofs, sub_gdofs, 1, 0, Comm);
  delete [] sub_gdofs;
  /*
  char fname[101]; int NumIndices, *Indices;
  sprintf(fname,"%s%d","test", MyPID);
  sprintf(fname, "%s.dat", fname);
  ofstream ffout;
  ffout.open(fname);
  ffout << ASt_red->NumMyRows() << endl;
  for (i=0; i<ASt_red->NumMyRows(); i++) ffout << ASt_red->GRID(i) << endl;
  ffout << ASt_red->NumMyCols() << endl;
  for (i=0; i<ASt_red->NumMyCols(); i++) ffout << ASt_red->GCID(i) << endl;
  ffout << RowMap->NumMyElements() << endl;
  for (i=0; i<RowMap->NumMyElements(); i++) ffout << RowMap->GID(i) << endl;
  ffout << ASt_red->Graph().MaxNumIndices() << endl;
  for (i=0; i<ASt_red->Graph().NumMyRows(); i++) {
    ASt_red->Graph().ExtractMyRowView(i, NumIndices, Indices);
    ffout << NumIndices << endl;
    for (int j=0; j<NumIndices; j++) ffout << Indices[j] << endl;
  }
  ffout.close();
  */
  for (i=0; i<= overlap; i++) {
    Importer = new Epetra_Import(*RowMap, ASt_red->RowMap());
    OverlapGraph = new Epetra_CrsGraph(Copy, *RowMap, 0);
    OverlapGraph->Import(ASt_red->Graph(), *Importer, Insert);
    OverlapGraph->FillComplete(ASt_red->DomainMap(), ASt_red->RangeMap());
    if (i < overlap) {
      delete RowMap;
      RowMap = new Epetra_BlockMap(OverlapGraph->ColMap());
      delete OverlapGraph;
      delete Importer;
    }
  }
  delete RowMap;
  Epetra_BlockMap rows = OverlapGraph->RowMap();
  delete OverlapGraph;
  delete Importer;
  int nrow = rows.NumMyElements();
  int *MyGlobalElements = new int[nrow];
  rows.MyGlobalElements(MyGlobalElements);
  Epetra_Map RowMap2(-1, nrow, MyGlobalElements, 0, Comm);
  delete [] MyGlobalElements;
  Epetra_Import Importer2(RowMap2, ASt_red->RowMap());
  Epetra_Import Importer3(RowMap2, AStandard->RowMap());
  AOverlap = new Epetra_CrsMatrix(Copy, RowMap2, RowMap2, 0);
  LDOverlap = new Epetra_IntVector(RowMap2);
  if (GNStandard) GNOverlap = new Epetra_IntVector(RowMap2);
  COverlap = new Epetra_MultiVector(RowMap2, 3);
  AOverlap->Import(*ASt_red, Importer2, Insert);
  AOverlap->FillComplete(ASt_red->DomainMap(), ASt_red->RangeMap());

  /*
  char fname[101];
  sprintf(fname,"%s%d","ASt_red", MyPID);
  sprintf(fname, "%s.dat", fname);
  spmat_datfile(*ASt_red, fname, 2);
  */

  LDOverlap->Import(*LDStandard, Importer3, Insert);
  if (GNStandard) GNOverlap->Import(*GNStandard, Importer3, Insert);
  COverlap->Import(*CStandard, Importer3, Insert);
  assert (COverlap->ConstantStride() == true);
  //  cout << *AOverlap << endl;
  // symmetry check
  //
  /*
  int j, k, dof, NumEntries, *Indices;
  double max_val(0), *Values, val1(0), val2(0);
  int N = AOverlap->NumMyRows();
  int *count_loc = new int[N]; myzero(count_loc, N);
  for (i=0; i<N; i++) {
    AOverlap->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) count_loc[Indices[j]]++;
  }
  int *rowbeg_loc = new int[N+1]; rowbeg_loc[0] = 0;
  for (i=0; i<N; i++) rowbeg_loc[i+1] = rowbeg_loc[i] + count_loc[i];
  int *colidx_loc = new int[rowbeg_loc[N]];
  double *val_loc = new double[rowbeg_loc[N]];
  myzero(count_loc, N);
  for (i=0; i<N; i++) {
    AOverlap->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      dof = Indices[j];
      colidx_loc[rowbeg_loc[dof] + count_loc[dof]] = i;
      val_loc[   rowbeg_loc[dof] + count_loc[dof]] = Values[j];
      count_loc[dof]++;
    }
  }
  delete [] count_loc;
  for (i=0; i<N; i++) {
    AOverlap->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      for (k=rowbeg_loc[i]; k<rowbeg_loc[i+1]; k++) {
	if (colidx_loc[k] == Indices[j]) {
	  double delta = fabs(val_loc[k] - Values[j]);
	  if (delta > max_val) {
	    max_val = delta;
	    val1 = Values[j];
	    val2 = val_loc[k];
	  }
	  break;
	}
      }
    }
  }
  delete [] rowbeg_loc; delete [] colidx_loc; delete [] val_loc;
  double normA = AOverlap->NormInf();
  if (MyPID == 0) cout << "infinity norm of AOverlap = " << normA << endl;
  cout << "maximum symmetry error    = " << max_val << endl;
  cout << "values = " << val1 << " " << val2 << endl;
  */
}

void CLOP_solver::construct_subdomains()
{
  if (print_flag > 9) fout << "in construct_subdomains" << endl;
  int i, j, partition_option(0);
  ndof_overlap = AOverlap->NumMyRows();
  count1 = new int[ndof_overlap];
  imap = new int[ndof_overlap];

  if ((GNOverlap) && (partition_option == 0)) {
    CLOP_graph CG(AOverlap, GNOverlap, overlap, partition_option,     
		  atype, ndim, count1, imap, dofpart1, dofpart2, npart);
    /*
    double *coords;
    COverlap->ExtractView(&coords, &ndof_overlap);
    double *x, *y, *z;
    x = &coords[0]; y = &coords[ndof_overlap]; z = &coords[2*ndof_overlap];
    if (MyPID == 0) {
      for (i=0; i<npart; i++) {
	cout << "MyPID = " << MyPID 
	     << ", dof coordinates for overlapping subdomain " << i << endl;
	for (j=dofpart2[i]; j<dofpart2[i+1]; j++) {
	  int dof = dofpart1[j];
	  cout << x[dof] << " " << y[dof] << endl;
	}
      }
    }
    */
  }
  else {
    npart = 1;
    dofpart2 = new int[npart+1]; dofpart2[0] = 0;
    dofpart2[1] = ndof_overlap;
    dofpart1 = new int[dofpart2[npart]];
    for (i=0; i<ndof_overlap; i++) dofpart1[i] = i;
  }
}

int CLOP_solver::initialize_subdomains()
{
  if (print_flag > 9) fout << "in initialize_subdomains" << endl;
  int i, j, max_nnz(0), nnz, dof, gdof, ipres, ipres_max, *locdof;
  max_ndof = 0; gmres_flag = 0;
  //
  // determine if pressure degrees of freedom present
  //
  ipres = 0;
  if (LDStandard != 0) {
    LDStandard->ExtractView(&locdof);
    for (i=0; i<LDStandard->MyLength(); i++)
      if (locdof[i] == 7) ipres = 1;
  }
  Comm.MaxAll(&ipres, &ipres_max, 1);
  if (ipres_max > 0) {
    gmres_flag = 1;
    if (print_flag > 0) fout << "Note: pressure dofs present" << endl;
  }
  //
  // initialize variables to be used later
  //
  memset(imap, -1, ndof_overlap*sizeof(int));
  csdima = new int[npart];
  ExporterO2ST = new Epetra_Export(AOverlap->RowMap(), ASt_red->RowMap());
  ImporterST2O = new Epetra_Import(AOverlap->RowMap(), ASt_red->RowMap());
  //
  // determine maximum number of dofs and maximum number of nonzeros
  // in stiffness matrix of each subdomain and allocate working memory
  //
  Asub = new CLOP_sub[npart];
  for (i=0; i<npart; i++) {
    int ndof_part = dofpart2[i+1] - dofpart2[i];
    if (ndof_part > max_ndof) max_ndof = ndof_part;
    Asub[i].getmatrix_nnz(&dofpart1[dofpart2[i]], ndof_part, AOverlap, 
			  imap, &Comm, nnz);
    if (nnz > max_nnz) max_nnz = nnz;
  }
  int *rowbeg_work = new int[max_ndof+1];
  int *colidx_work = new int[max_nnz];
  double *A_work = new double[max_nnz]; myzero(A_work, max_nnz);
  //
  // determine starting subdomain number for each processor
  //
  Comm.ScanSum(&npart, &gpart0, 1); gpart0 -= npart;
  //
  // determine subdomains (both on and off-processor) associated with
  // each dof and store in Epetra Crs_Graph Overlap_Subs
  // (row i of Overlap_Subs contains global subdomain numbers containing
  // local overlap dof i)
  //
  Epetra_CrsGraph *Overlap_Subs;
  construct_Overlap_Subs(Overlap_Subs);
  //
  // flag dofs on subdomain boundaries and construct nsubdof array
  // on_sub_bound[i] = 0 if local dof i not on subdomain boundary
  //                 = 1 if local dof i is  on subdomain boundary
  // nsubdof[i]      = number of subdomains containing local dof i
  //
  unsigned char *on_sub_bound = new unsigned char[ndof_overlap];
  unsigned char *nsubdof      = new unsigned char[ndof_overlap];
  flag_sub_bound(on_sub_bound, nsubdof, Overlap_Subs);
  delete Overlap_Subs;
  /*
  if (MyPID == 1) {
    double *coords;
    COverlap->ExtractView(&coords, &ndof_overlap);
    double *x, *y, *z;
    x = &coords[0]; y = &coords[ndof_overlap]; z = &coords[2*ndof_overlap];
    for (i=0; i<ndof_overlap; i++) {
      cout << x[i] << " " << y[i] << " " << z[i] << " " <<
	1*nsubdof[i] << " " << 1*on_sub_bound[i] << endl;
    }
  }
  */
  //
  // code for exact solver follows
  //
  if (local_solver == 1) {
    int max_nrhs;
    if (atype == 1) max_nrhs = ndim + 1;
    if (atype == 2) {
      if (ndim == 2) max_nrhs = 6;
      if (ndim == 3) max_nrhs = 12;
    }
    if (atype == 3) max_nrhs = 9;
    //
    // factor stiffness matrix of each subdomain
    //
    if (print_flag > 9) fout << "factoring subdomain matrices" << endl;
    for (i=0; i<npart; i++) {
      Asub[i].factormatrix(AOverlap, imap, rowbeg_work, colidx_work, A_work);
    }
    //
    // construct partition of unity
    //
    int max_dim1 = max_ndof;
    if (max_nrhs > max_dim1) max_dim1 = max_nrhs;
    rhs_work = new double[max_dim1*max_nrhs];
    sol_work = new double[max_dim1*max_nrhs];
    tmp_work = new double[max_dim1*max_nrhs];
    int LWORK(250);
    double *WORK = new double[LWORK];
    double *Edof = new double[ndof_overlap]; myzero(Edof, ndof_overlap);
    int nneg, nneg_max, nneg_max_proc(0);
    for (i=0; i<npart; i++) {
      Asub[i].genpu(LDOverlap, COverlap, rhs_work, sol_work, tmp_work, 
		    atype, ndim, WORK, LWORK, Edof, nneg);
      if (nneg > nneg_max_proc) nneg_max_proc = nneg;
    }
    Comm.MaxAll(&nneg_max_proc, &nneg_max, 1);
    //  assert (nneg_max == 0);
    if (nneg_max > 0) {
      gmres_flag = 1;
      if (print_flag >= 0) {
	fout << "Warning: stiffness matrix is not positive definite " << endl;
	fout << "         gmres_flag set to 1" << endl;
      }
    }
    Epetra_Vector Edof_Overlap( View, AOverlap->RowMap(), Edof);
    //    cout << Edof_Overlap << endl;
    Epetra_Vector Edof_Standard(ASt_red->RowMap());
    //    cout << Edof_Overlap << endl;
    Edof_Standard.Export(Edof_Overlap, *ExporterO2ST, Add);
    //    cout << Edof_Standard << endl;
    Edof_Overlap.Import(Edof_Standard, *ImporterST2O, Insert);
    //    cout << Edof_Overlap << endl;
    for (i=0; i<npart; i++) Asub[i].normalpu(Edof, nsubdof);
    //
    // construct coarse space (pass 1)
    //
    ncdof_proc = 0; max_csdim = 0;
    for (i=0; i<npart; i++) {
      Asub[i].construct_coarse1(AOverlap, rhs_work, sol_work, tmp_work,
				rowbeg_work, colidx_work, A_work, imap, 
				nsubdof, csdima[i], ndof_rot);
      ncdof_proc += csdima[i];
      if (csdima[i] > max_csdim) max_csdim = csdima[i];
    }
    //    cout << "MyPID, ncdof_proc = " << MyPID << " " << ncdof_proc << endl;
    cs_local = new int[ncdof_proc];
    xcent = new double[npart]; 
    ycent = new double[npart]; 
    zcent = new double[npart];
    ncdof_proc = 0;
    for (i=0; i<npart; i++) {
      Asub[i].get_cdof(&cs_local[ncdof_proc], csdima[i], xcent[i], ycent[i], 
		       zcent[i]);
      ncdof_proc += csdima[i];
    }
    /*
    cout << "subdomain centroids" << endl;
    for (i=0; i<npart; i++) cout << xcent[i] << " " << ycent[i] << " " 
				 << zcent[i] << endl;
    cout << "cs_local = "; for (i=0; i<ncdof_proc; i++) cout << cs_local[i] 
	 << " ";
    cout << endl;
    */
    //
    // correct coarse space if rotational dofs present in model
    //
    int max_ndof_rot;
    Comm.MaxAll(&ndof_rot, &max_ndof_rot, 1);
    //    cout << "max_ndof_rot = " << max_ndof_rot << endl;
    if (max_ndof_rot > 0) {
      if ((atype == 2) && (ndim == 3)) {
	for (int rbm=1; rbm<7; rbm++) {
	  correct_shape_ela(rbm, Edof, Edof_Overlap, Edof_Standard, nsubdof);
	}
      }
      if (atype == 3) {
	for (int rbm=1; rbm<4; rbm++) {
	  correct_shape_dkt(rbm, Edof, Edof_Overlap, Edof_Standard, nsubdof);
	}
      }
    }
    //
    // perform static condensation energy minimization for coarse space
    //
    for (i=0; i<npart; i++) {
      Asub[i].statcond(nsubdof, on_sub_bound, imap, rowbeg_work, colidx_work, 
		       A_work, rhs_work, sol_work, tmp_work, AOverlap);
    }
    delete [] imap;
    delete [] nsubdof;
    delete [] on_sub_bound;
    delete [] rowbeg_work;
    delete [] colidx_work;
    delete [] A_work;
    delete [] rhs_work; rhs_work = 0;
    delete [] sol_work; sol_work = 0;
    delete [] tmp_work; tmp_work = 0;
    delete [] WORK;
    delete [] Edof;
    //
    // assemble coarse interpolation matrix
    //
    assemble_Phi();
    return(0);
  }
  return(0);
}

void CLOP_solver::construct_Overlap_Subs(Epetra_CrsGraph* & Overlap_Subs)
{
  int i, j, dof;
  myzero(count1, ndof_overlap);
  for (i=0; i<npart; i++) {
    for (j=dofpart2[i]; j<dofpart2[i+1]; j++) count1[dofpart1[j]]++;
  }
  int *sub2 = new int[ndof_overlap+1]; sub2[0] = 0;
  for (i=0; i<ndof_overlap; i++) sub2[i+1] = sub2[i] + count1[i];
  int *sub1 = new int[sub2[ndof_overlap]];
  myzero(count1, ndof_overlap);
  int *global_cols = new int[npart];
  for (i=0; i<npart; i++) {
    global_cols[i] = gpart0 + i;
    for (j=dofpart2[i]; j<dofpart2[i+1]; j++) {
      dof = dofpart1[j];
      sub1[sub2[dof] + count1[dof]] = i;
      count1[dof]++;
    }
  }
  Epetra_Map ColMap(-1, npart, global_cols, 0, Comm);
  Epetra_CrsGraph Overlap_Graph(Copy, AOverlap->RowMap(), ColMap, count1);
  for (i=0; i<ndof_overlap; i++) {
    Overlap_Graph.InsertMyIndices(i, count1[i], &sub1[sub2[i]]);
  }
  delete [] sub1; 
  delete [] sub2; 
  delete [] global_cols;
  Overlap_Graph.FillComplete(ColMap, ASt_red->RowMap());
  Epetra_CrsGraph Standard_Graph(Copy, ASt_red->RowMap(), 0);
  /*
  char fname[101];
  sprintf(fname,"%s%d","test_g", MyPID);
  sprintf(fname, "%s.dat", fname);
  ofstream ffout;
  ffout.open(fname);
  int NumEntries, *Indices;
  ffout << Overlap_Graph.MaxNumIndices() << endl;
  for (i=0; i<Overlap_Graph.NumMyRows(); i++) {
    Overlap_Graph.ExtractMyRowView(i, NumEntries, Indices);
    ffout << NumEntries << endl;
    for (j=0; j<NumEntries; j++) ffout << Indices[j] << endl;
  }
  ffout.close();
  */
  Standard_Graph.Export(Overlap_Graph, *ExporterO2ST, Insert);
  Standard_Graph.FillComplete(ColMap, ASt_red->RowMap());
  Overlap_Subs = new Epetra_CrsGraph(Copy, AOverlap->RowMap(), 0);
  Overlap_Subs->Import(Standard_Graph, *ImporterST2O, Insert);
  Overlap_Subs->FillComplete(ColMap, ASt_red->RowMap());
}

void CLOP_solver::flag_sub_bound(unsigned char on_sub_bound[], 
   unsigned char nsubdof[], Epetra_CrsGraph* & Overlap_Subs)
{
  int i, j, k, m, dof, NumEntries, *Indices, iflag, isub;
  int NumIndices1, NumIndices2, *Indices1, *Indices2;
  double *Values;
  myzero(on_sub_bound, ndof_overlap);
  for (i=0; i<ndof_overlap; i++) {
    Overlap_Subs->ExtractMyRowView(i, NumIndices1, Indices1);
    nsubdof[i] = NumIndices1;
    AOverlap->ExtractMyRowView(i, NumEntries, Values, Indices);
    assert(NumEntries > 0);
    for (j=0; j<NumEntries; j++) {
      dof = Indices[j];
      Overlap_Subs->ExtractMyRowView(dof, NumIndices2, Indices2);
      for (k=0; k<NumIndices2; k++) {
	isub = Indices2[k];
	iflag = 1;
	for (m=0; m<NumIndices1; m++) {
	  if (Indices1[m] == isub) {
	    iflag = 0;
	    break;
	  }
	}
	if (iflag == 1) break;
      }
      if (iflag == 1) break;
    }
    if (iflag == 1) on_sub_bound[i] = 1;
  }
  //
  // sum contributions to on_sub_bound over all processors
  //
  for (i=0; i<ndof_overlap; i++) count1[i] = on_sub_bound[i];
  int nSt_red = ASt_red->NumMyRows();
  int *ivec_St_red = new int[nSt_red]; myzero(ivec_St_red, nSt_red);
  Epetra_IntVector osb_St_red(View, ASt_red->RowMap(), ivec_St_red);
  Epetra_IntVector osb_Overlap(View, AOverlap->RowMap(), count1);
  osb_St_red.Export(osb_Overlap, *ExporterO2ST, Add);
  osb_Overlap.Import(osb_St_red, *ImporterST2O, Insert);
  for (i=0; i<ndof_overlap; i++) if (count1[i] > 0) on_sub_bound[i] = 1;
  delete [] ivec_St_red;
}

void CLOP_solver::correct_shape_ela(int rbm, double Edof[],
	  Epetra_Vector & Edof_Overlap, Epetra_Vector & Edof_Standard, 
	  unsigned char nsubdof[])
{
  int i, *locdof;
  LDOverlap->ExtractView(&locdof);
  myzero(Edof, ndof_overlap);
  for (i=0; i<npart; i++) Asub[i].sum_scalar_multiply(Edof, rbm, 1.0);
  if (rbm > 3) {
    if (rbm == 4) {
      for (i=0; i<npart; i++) {
	Asub[i].sum_scalar_multiply(Edof, 2, -zcent[i]);
	Asub[i].sum_scalar_multiply(Edof, 3,  ycent[i]);
      }
    }
    if (rbm == 5) {
      for (i=0; i<npart; i++) {
	Asub[i].sum_scalar_multiply(Edof, 3, -xcent[i]);
	Asub[i].sum_scalar_multiply(Edof, 1,  zcent[i]);
      }
    }
    if (rbm == 6) {
      for (i=0; i<npart; i++) {
	Asub[i].sum_scalar_multiply(Edof, 1, -ycent[i]);
	Asub[i].sum_scalar_multiply(Edof, 2,  xcent[i]);
      }
    }
  }
  Edof_Standard.PutScalar(0.0);
  Edof_Standard.Export(Edof_Overlap, *ExporterO2ST, Add);
  Edof_Overlap.Import(Edof_Standard, *ImporterST2O, Insert);
  for (i=0; i<npart; i++) Asub[i].correct(Edof, rbm, nsubdof, locdof, 3, 7);
}

void CLOP_solver::correct_shape_dkt(int rbm, double Edof[],
       Epetra_Vector & Edof_Overlap, Epetra_Vector & Edof_Standard, 
       unsigned char nsubdof[])
{
  int i, *locdof;
  LDOverlap->ExtractView(&locdof);
  myzero(Edof, ndof_overlap);
  for (i=0; i<npart; i++) Asub[i].sum_scalar_multiply(Edof, rbm, 1.0);
  if (rbm > 1) {
    if (rbm == 2) {
      for (i=0; i<npart; i++) {
	Asub[i].sum_scalar_multiply(Edof, 1,  ycent[i]);
      }
    }
    if (rbm == 3) {
      for (i=0; i<npart; i++) {
	Asub[i].sum_scalar_multiply(Edof, 1, -xcent[i]);
      }
    }
  }
  Edof_Standard.PutScalar(0.0);
  Edof_Standard.Export(Edof_Overlap, *ExporterO2ST, Add);
  Edof_Overlap.Import(Edof_Standard, *ImporterST2O, Insert);
  for (i=0; i<npart; i++) Asub[i].correct(Edof, rbm, nsubdof, locdof, 1, 4);
}

void CLOP_solver::assemble_Phi()
{
  if (print_flag > 9) fout << "in assemble_Phi" << endl;
  int i, j, k, ldof, row, ibeg, jbeg, col, ierr;
  //
  // determine starting coarse dof number for each processor
  //
  Comm.ScanSum(&ncdof_proc, &gcdof0, 1); gcdof0 -= ncdof_proc;
  cdof_proc = new int[ncdof_proc];
  for (i=0; i<ncdof_proc; i++) cdof_proc[i] = gcdof0 + i;
  int *NumEntriesPerRow = NULL;
  row = 0;
  NumEntriesPerRow = new int[ncdof_proc];
  for (i=0; i<npart; i++) {
    for (j=0; j<csdima[i]; j++) {
      NumEntriesPerRow[row] = dofpart2[i+1] - dofpart2[i];
      row++;
    }
  }
  //
  // construct Epetra CrsMatrix for transpose of Phi (PhiT)
  //
  Epetra_Map PhiTRowMap(-1, ncdof_proc, cdof_proc, 0, Comm);
  PhiT = new Epetra_CrsMatrix(Copy, PhiTRowMap, AOverlap->RowMap(),
			      NumEntriesPerRow); 
  double *Phi_ptr;
  row = 0;
  for (i=0; i<npart; i++) {
    Asub[i].get_Phi_ptr(Phi_ptr);
    ibeg = 0;
    for (j=0; j<csdima[i]; j++) {
      ierr = PhiT->InsertMyValues(row, NumEntriesPerRow[row], &Phi_ptr[ibeg],
				  &dofpart1[dofpart2[i]]);
      assert (ierr == 0);
      ibeg += NumEntriesPerRow[row];
      row++;
    }
  }
  PhiT->FillComplete(ASt_red->RowMap(), PhiTRowMap);
  //  cout << *PhiT << endl;
  delete [] NumEntriesPerRow;
  //
  // construct Epetra CrsMatrix for Phi (Phi_overlap, on-processor coarse 
  // dofs only, this is the transpose of PhiT)
  //
  NumEntriesPerRow = new int[ndof_overlap]; 
  myzero(NumEntriesPerRow, ndof_overlap);
  for (i=0; i<npart; i++) {
    for (j=dofpart2[i]; j<dofpart2[i+1]; j++) {
      ldof = dofpart1[j];
      NumEntriesPerRow[ldof] += csdima[i];
    }
  }
  Epetra_CrsMatrix Phi_Overlap(Copy, AOverlap->RowMap(), PhiTRowMap,
			       NumEntriesPerRow);
  int *rowbeg_Phi = new int[ndof_overlap+1]; rowbeg_Phi[0] = 0;
  for (i=0; i<ndof_overlap; i++) {
    rowbeg_Phi[i+1] = rowbeg_Phi[i] + NumEntriesPerRow[i];
  }
  int nnz_Phi = rowbeg_Phi[ndof_overlap];
  int *colidx_Phi = new int[nnz_Phi];
  double *phivec = new double[nnz_Phi];
  myzero(count1, ndof_overlap);
  col = 0;
  for (i=0; i<npart; i++) {
    Asub[i].get_Phi_ptr(Phi_ptr);
    jbeg = 0;
    for (j=0; j<csdima[i]; j++) {
      for (k=dofpart2[i]; k<dofpart2[i+1]; k++) {
	ldof = dofpart1[k];
	colidx_Phi[rowbeg_Phi[ldof] + count1[ldof]] = col;
	phivec[rowbeg_Phi[ldof] + count1[ldof]] = Phi_ptr[jbeg];
	count1[ldof]++;
	jbeg++;
      }
      col++;
    }
  }
  for (i=0; i<ndof_overlap; i++) {
    ierr = Phi_Overlap.InsertMyValues(i, NumEntriesPerRow[i], 
        &phivec[rowbeg_Phi[i]], &colidx_Phi[rowbeg_Phi[i]]);
    assert (ierr == 0);
  }
  ierr = Phi_Overlap.FillComplete(PhiTRowMap, ASt_red->RowMap());
  assert (ierr == 0);
  //  cout << Phi_Overlap << endl;
  Phi = new Epetra_CrsMatrix(Copy, ASt_red->RowMap(), 0);
  //  cout << Phi_Overlap.RowMap() << endl;
  //  cout << Phi->RowMap() << endl;
  //  cout << *ExporterO2ST << endl;
  Epetra_Export newExporter(Phi_Overlap.RowMap(), Phi->RowMap());

  //  assert (newExporter.SourceMap().SameAs(ExporterO2ST->SourceMap()));
  //  assert (newExporter.TargetMap().SameAs(ExporterO2ST->TargetMap()));
  //  Comm.Barrier();
  /*
  {
  char fname[101];
  sprintf(fname,"%s%d","test", MyPID);
  sprintf(fname, "%s.dat", fname);
  ofstream ffout;
  ffout.open(fname);
  int N, NumEntries, *Indices;
  double *Values;
  N = ExporterO2ST->SourceMap().NumMyElements();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << ExporterO2ST->SourceMap().GID(i) << endl;

  N = ExporterO2ST->TargetMap().NumMyElements();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << ExporterO2ST->TargetMap().GID(i) << endl;

  N = PhiTRowMap.NumMyElements();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << PhiTRowMap.GID(i) << endl;

  ffout << Phi_Overlap.MaxNumEntries() << endl;
  for (i=0; i<Phi_Overlap.NumMyRows(); i++) {
    Phi_Overlap.ExtractMyRowView(i, NumEntries, Values, Indices);
    ffout << NumEntries << endl;
    for (j=0; j<NumEntries; j++) ffout << Indices[j] << " " 
				      << Values[j] << endl;
  }
  ffout.close();
  }
  */
  Phi->Export(Phi_Overlap, *ExporterO2ST, Insert);
  //  Phi->Export(Phi_Overlap, newExporter, Insert);
  Phi->FillComplete(PhiTRowMap, ASt_red->RowMap());
  int *Indices, NumEntries;
  double *Values;
  /*
  for (i=0; i<Phi->NumMyRows(); i++) {
    double dsum = 0;
    Phi->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) dsum += Values[j];
    if (MyPID == 3) cout << Phi->GRID(i) << ": " << dsum << endl;
  }
  */
  //  cout << *Phi << endl;
  //  cout << *Phi << endl;
  delete [] NumEntriesPerRow;
  delete [] rowbeg_Phi;
  delete [] colidx_Phi;
  delete [] phivec;
  delete [] count1;
}

void CLOP_solver::calculate_coarse_stiff()
{
  if (print_flag > 9) fout << "in calculate_coarse_stiff" << endl;
  Epetra_CrsMatrix *KPhi;
  EpetraExtCD::MatrixMatrix::Multiply(*ASt_red, *Phi, KPhi);
  EpetraExtCD::MatrixMatrix::Multiply(*PhiT, *KPhi, Kc);
  delete KPhi;
}

void CLOP_solver::gather_coarse_stiff()
{
  if (print_flag > 9) fout << "in gather_coarse_stiff" << endl;
  int i, ncdof_solver;
  assert (Kc->NumGlobalCols() == Kc->NumGlobalRows());
  ncdof = Kc->NumGlobalRows();
  if (MyPID == 0) ncdof_solver = ncdof;
  if (MyPID != 0) ncdof_solver = 0;
  int *allrows = NULL;
  allrows = new int[ncdof_solver];
  for (i=0; i<ncdof_solver; i++) allrows[i] = i;
  RowMap_coarse = new Epetra_Map(ncdof, ncdof_solver, allrows, 0, Comm);
  delete [] allrows;
  Kc_gathered = new Epetra_CrsMatrix(Copy, *RowMap_coarse, 0);
  Importer_coarse = new Epetra_Import(*RowMap_coarse, Kc->RowMap());
  Kc_gathered->Import(*Kc, *Importer_coarse, Insert);
  Kc_gathered->FillComplete();
  //  cout << *Kc_gathered << endl;
  /*
  char fname[101];
  sprintf(fname,"%s","coarse_mat.dat");
  spmat_datfile(*Kc_gathered, fname, 1);
  */
}

void CLOP_solver::factor_coarse_stiff()
{
  if (print_flag > 9) fout << "factoring coarse matrix" << endl;
  int i, j;
  Kc_gathered->MakeDataContiguous();
  int *rowbeg_KC, *colidx_KC, NumEntries;
  double *KC;
  int nnz_KC = Kc_gathered->NumMyNonzeros();
  if (MyPID == 0) {
    rowbeg_KC = new int[ncdof+1]; rowbeg_KC[0] = 0;
    for (i=0; i<ncdof; i++) {
      Kc_gathered->ExtractMyRowView(i, NumEntries, KC, colidx_KC);
      rowbeg_KC[i+1] = rowbeg_KC[i] + NumEntries;
    }
    if (ncdof > 0) Kc_gathered->ExtractMyRowView(0, NumEntries, KC, colidx_KC);
    /*
    cout << "ncdof, nnz_KC = " << ncdof << " " << nnz_KC << endl;
    for (i=0; i<ncdof; i++) {
      for (j=rowbeg_KC[i]; j<rowbeg_KC[i+1]; j++) {
	cout << i << " " << colidx_KC[j] << " " << KC[j] << endl;
      }
    }
    ofstream ffout;
    ffout.open("coarse_rcv.dat");
    ffout << ncdof << endl;
    ffout << nnz_KC << endl;
    for (i=0; i<=ncdof; i++) ffout << rowbeg_KC[i] << endl;
    for (i=0; i<nnz_KC; i++) ffout << colidx_KC[i] << " " << KC[i] << endl;
    ffout.close();
    */

    Kc_fac = new sparse_lu();
    if (ncdof <= ndof_global_red)
      Kc_fac->factor(ncdof, nnz_KC, rowbeg_KC, colidx_KC, KC);
    delete [] rowbeg_KC;
  }
}

void CLOP_solver::solve_init()
{
  int ncdof_sol = Kc_gathered->NumMyRows();
  rhs_coarse  = new double[ncdof_sol];
  sol_coarse  = new double[ncdof_sol];
  temp_coarse = new double[ncdof_sol];
  r_overlap = new double[ndof_overlap];
  z_overlap = new double[ndof_overlap];
  rhs_work = new double[max_ndof];
  sol_work = new double[max_ndof];
  tmp_work = new double[max_ndof];
  rcurra = new double[maxiter+2];
  rhoa = new double[maxiter+2];
  betaa = new double[maxiter+2];
  pApa = new double[maxiter+2];
  Dtri = new double[maxiter+2];
  Etri = new double[maxiter+2];
  econa = new double[maxiter+2];
  rSt_red = new Epetra_Vector(ASt_red->RowMap());
  zSt_red = new Epetra_Vector(ASt_red->RowMap());
  pSt_red = new Epetra_Vector(ASt_red->RowMap());
  vSt_red = new Epetra_Vector(ASt_red->RowMap());
  if (ncon_global == 0) {
    uSt_red = new Epetra_Vector(ASt_red->RowMap());
    vStand = new Epetra_Vector(AStandard->RowMap());
  }
  ApSt_red = new Epetra_Vector(ASt_red->RowMap());
  wStand = new Epetra_Vector(AStandard->RowMap());
  rOverlap = new Epetra_Vector(View, AOverlap->RowMap(), r_overlap);
  zOverlap = new Epetra_Vector(View, AOverlap->RowMap(), z_overlap);
  PhiTr  = new Epetra_Vector(PhiT->RowMap());
  PhiTr_gathered = new Epetra_Vector(View, Kc_gathered->RowMap(), rhs_coarse);
  CSol_gathered = new Epetra_Vector(View, Kc_gathered->RowMap(), sol_coarse);
  int nvec = max_orthog; if (nvec == 0) nvec = 1;
  if (gmres_flag == 0) {
    AP_matrix = new Epetra_MultiVector(ASt_red->RowMap(), nvec);
    P_matrix  = new Epetra_MultiVector(ASt_red->RowMap(), nvec);
    gSt_red = new Epetra_Vector(ASt_red->RowMap());
    Epetra_LocalMap LocalMap(max_orthog, 0, Comm);
    pAp_vec   = new double[nvec];
    ortho_vec = new double[nvec];
    ortho_sum = new double[nvec];
    PAP = new double[4*nvec];
    PAP_sum = new double[nvec*nvec];
    PAP_store = new double[nvec*nvec];
    IPIV = new int[nvec];
    n_orthog_used = 0;
  }
  if (gmres_flag == 1) {
    int nn = ASt_red->NumMyRows();
    int m = nn*(maxiter+1);
    VV = new double[m]; myzero(VV, m);
    m = maxiter*maxiter;
    RR = new double[m]; myzero(RR, m);
    m = maxiter+1;
    HH = new double[m]; myzero(HH, m);
    zz = new double[m]; myzero(zz, m);
    cc = new double[m]; myzero(cc, m);
    ss = new double[m]; myzero(ss, m);
    norms = new double[m]; myzero(norms, m);
    gmres_vec = new double[m];
    gmres_sum = new double[m];
  }
}

void CLOP_solver::solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
			int & num_iter, int & solver_status)
{
  int pcg_status, gmres_status;
  double starttime, endtime;
  Comm.Barrier();
  if (MyPID == 0) starttime = MPI_Wtime();
  if (print_flag >= 0) fout.open("CLOP_solver.data", ios::app);
  if (gmres_flag == 0) {
    pcg_solve(uStand, fStand, num_iter, pcg_status);
    solver_status = pcg_status;
  }
  if (gmres_flag == 1) {
    gmres_solve(uStand, fStand, num_iter, gmres_status);
    solver_status = gmres_status;
  }
  Comm.Barrier();
  if (print_flag > 0) {
    endtime = MPI_Wtime();
    fout << "elapsed time for clop solver solve = " << endtime-starttime
	 << " seconds" << endl;
  }
  if (print_flag >= 0) fout.close();
}

void CLOP_solver::pcg_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
			    int & num_iter, int & pcg_status)
{
  int i, j, iflag(1), pcg_iter, max_iter_new;
  double dprod, rorig, beta, roldzold(0), pAp, alpha, rcurr, ractual;
  double norm_rconstraint, norm_conerror;
  orthog_option = 2;
  //
  // determine reduced residual
  //
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *rSt_red);
  if (ncon_global == 0) rSt_red->Update(1.0, *fStand, 0.0);
  //
  // calculate initial residual
  //
  rSt_red->Dot(*rSt_red, &dprod);
  rorig = sqrt(dprod);
  if (print_flag > 0) fout << "rorig = " << rorig << endl;
  rcurra[0] = rorig;

  //
  // initial coarse grid correction and residual update if multiplicative
  // coarse grid correction selected
  //
  pre_type_coarse = 1; // multiplicative coarse grid correction
  if (pre_type_coarse == 1) {
    coarse_correction(rSt_red, uSt_red);
    mat_vec_prod(uSt_red, ApSt_red);
    rSt_red->Update(-1.0, *ApSt_red, 1.0);
  }
  //
  // calculate correction from stored search directions and update
  // residual
  //
  pcg_status = 1;
  max_iter_new = maxiter;
  if (n_orthog_used > 0) {
    search_correction(rSt_red, pSt_red, 0.0);
    mat_vec_prod(pSt_red, ApSt_red);
    rSt_red->Update(-1.0, *ApSt_red, 1.0);
    rSt_red->Dot(*rSt_red, &dprod);
    if (MyPID == 0) {
      cout << "residual after initial orthog correction = " << sqrt(dprod) 
	   << endl << "number of search directions used = "
	   << n_orthog_used << endl;
    }
    uSt_red->Update(1.0, *pSt_red, 1.0);
    if ((iflag > 0) && (sqrt(dprod)/rorig <= solver_tol)) {
      pcg_status = max_iter_new = num_iter = 0;
    }
  }
  //
  // pcg iterations
  //
  pre_type_orthog = 0; if (n_orthog_used == max_orthog) pre_type_orthog = 1;
  pcg_iter = 0;
  for (int iter=0; iter<max_iter_new; iter++) {
    //    if (MyPID == 0) cout << "iteration " << iter+1 << " of maxiter = "
    //			 << maxiter << endl;
    apply_preconditioner(rSt_red, zSt_red);

    if (pre_type_orthog == 0) { // project off stored search directions
      pSt_red->Update(1.0, *zSt_red, 0.0);
      if (n_orthog_used > 0) remove_projection_search(pSt_red);
      rSt_red->Dot(*pSt_red, &dprod);
    }
    else { // standard preconditioned conjugate gradient stuff 
      rSt_red->Dot(*zSt_red, &dprod);
      assert (dprod >= 0);
      rhoa[pcg_iter] = sqrt(fabs(dprod));
      if (pcg_iter == 0) {
	beta = 0;
	pSt_red->Update(1.0, *zSt_red, 0.0);
      }
      else {
	beta = dprod/roldzold;
	pSt_red->Update(1.0, *zSt_red, beta);
      }
      betaa[pcg_iter] = beta;
      roldzold = dprod;
      pcg_iter++;
    }
    //
    // determine change in residual caused by p
    //
    mat_vec_prod(pSt_red, ApSt_red);
    pSt_red->Dot(*ApSt_red, &pAp);
    pApa[iter] = pAp;
    alpha = dprod/pAp;
    //
    // store search direction and stiffness matrix times search direction
    //
    if (n_orthog_used == max_orthog) pre_type_orthog = 1;
    store_search(pAp);
    //
    // update uStand and rStand
    //
    uSt_red->Update( alpha,  *pSt_red, 1.0);
    rSt_red->Update(-alpha, *ApSt_red, 1.0);
    rSt_red->Dot(*rSt_red, &rcurr);
    rcurr = sqrt(rcurr);
    rcurra[iter+1] = rcurr;
    num_iter = iter+1;
    //    if (MyPID == 0) cout << "alpha, dprod, rtol = " << alpha << " " 
    //			 << dprod << " " << rcurr/rorig << endl;
    if ((iflag > 0) && (rcurr/rorig <= solver_tol)) {
      pcg_status = 0;
      break;
    }
  }
  if (ncon_global >  0) {
    Tran->Multiply(false, *uSt_red, *uStand);
    Tran->Multiply(true, *fStand, *rSt_red);
  }
  if (ncon_global == 0) {
    uStand->Update(1.0, *uSt_red, 0.0);
    rSt_red->Update(1.0, *fStand, 0.0);
  }
  //
  // calculate actual residual for reduced system
  //
  ASt_red->Multiply(false, *uSt_red, *zSt_red);
  zSt_red->Update(1.0, *rSt_red, -1.0);
  zSt_red->Dot(*zSt_red, &ractual);
  ractual = sqrt(ractual);
  //
  // calculate Lagrange multipliers and check residual for full system
  //
  if (ncon_global > 0) {
    AStandard->Multiply(false, *uStand, *vStand);
    vStand->Update(1.0, *fStand, -1.0);
    calculate_multipliers(uStand, norm_rconstraint, norm_conerror);
  }
  if (MyPID == 0) {
    //    cout << "rorig                 = " << rorig << endl;
    if (print_flag > 0) {
      if (num_iter > 0) fout << "rcurr(recursive)      = " << rcurr << endl;
      fout << "rcurr(actual)         = " << ractual << endl;
      if (ncon_global > 0) {
	fout << "rcurr(constraint)     = " << norm_rconstraint << endl;
	fout << "constraint error norm = " << norm_conerror << endl;
      }
      fout << "number of iterations  = " << num_iter << endl;
      fout << "solver tolerance      = " << solver_tol << endl;
    }
    if ((print_flag == 2) || (print_flag == 12)) {
      fout << "condition # estimate      relative residual" 
	   << "   iteration" << endl;
      if (pcg_iter == num_iter) calculate_condition(pcg_iter);
      fout << setiosflags(ios::scientific | ios::uppercase);
      for (i=0; i<num_iter; i++) {
	double ee = 0;
	if (pcg_iter == num_iter) ee = econa[i]; 
	fout << " " 
	     << setw(17) << setprecision(10) << ee 
	     << "       " 
	     << setw(17) << setprecision(10) << rcurra[i+1]/rorig
	     << "        " 
	     << i+1 << endl;
      }
    }
    fout << resetiosflags(ios::scientific);
    fout << resetiosflags(ios::uppercase);
    fout << setprecision(6);
  }
}

void CLOP_solver::gmres_solve(Epetra_Vector* uStand, 
          const Epetra_Vector* fStand, int & num_iter, int & gmres_status)
{
  int i, j, iflag(1), gmres_iter;
  double dprod, rorig, beta, normb, rcurr, ractual;
  double norm_rconstraint, norm_conerror;
  //
  // determine reduced residual
  //
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *rSt_red);
  if (ncon_global == 0) rSt_red->Update(1.0, *fStand, 0.0);
  //
  // calculate initial residual
  //
  rSt_red->Dot(*rSt_red, &dprod);
  rorig = sqrt(dprod);
  if (MyPID == 0) cout << "rorig = " << rorig << endl;
  uSt_red->PutScalar(0);
  //
  // gmres iterations
  //
  normb = rorig;
  gmres_status = 1;
  int n = ASt_red->NumMyRows();
  double *vals = rSt_red->Values();
  for (i=0; i<n; i++) {
    vals[i] /= normb;
    VV[i] = vals[i];
  }
  pre_type_orthog = 0; pre_type_coarse = 1;
  for (gmres_iter=0; gmres_iter<maxiter; gmres_iter++) {
    //    if (MyPID == 0) cout << "iteration " << gmres_iter+1 
    //			 << " of maxiter = " << maxiter << endl;
    apply_preconditioner(rSt_red, zSt_red);
    //
    // gmres stuff
    //
    mat_vec_prod(zSt_red, rSt_red);
    //
    // two steps of classical Gram-Schmidt
    //
    two_steps_CGS(gmres_iter, rSt_red);
    //
    // Hessengberg QR using two by two rotations
    //
    hessenberg_qr(gmres_iter);
    rcurr = normb*fabs(norms[gmres_iter]);
    num_iter = gmres_iter+1;
    if ((iflag > 0) && (rcurr/rorig <= solver_tol)) {
      gmres_status = 0;
      break;
    }
    memcpy(vals, &VV[n*num_iter], n*sizeof(double));
  }
  //
  // construct the solution
  //
  construct_solution(num_iter, normb);
  if (ncon_global >  0) {
    Tran->Multiply(false, *uSt_red, *uStand);
    Tran->Multiply(true, *fStand, *rSt_red);
  }
  if (ncon_global == 0) {
    uStand->Update(1.0, *uSt_red, 0.0);
    rSt_red->Update(1.0, *fStand, 0.0);
  }
  //
  // calculate actual residual for reduced system
  //
  ASt_red->Multiply(false, *uSt_red, *zSt_red);
  zSt_red->Update(1.0, *rSt_red, -1.0);
  zSt_red->Dot(*zSt_red, &ractual);
  ractual = sqrt(ractual);
  //
  // calculate Lagrange multipliers and check residual for full system
  //
  if (ncon_global > 0) {
    AStandard->Multiply(false, *uStand, *vStand);
    vStand->Update(1.0, *fStand, -1.0);
    calculate_multipliers(uStand, norm_rconstraint, norm_conerror);
  }
  if (MyPID == 0) {
    //    cout << "rorig                 = " << rorig << endl;
    cout << "rcurr(gmres)          = " << rcurr << endl;
    cout << "rcurr(actual)         = " << ractual << endl;
    if (ncon_global > 0) {
      cout << "rcurr(constraint)     = " << norm_rconstraint << endl;
      cout << "constraint error norm = " << norm_conerror << endl;
    }
    cout << "number of iterations  = " << num_iter << endl;
    cout << "solver tolerance      = " << solver_tol << endl;
    if (iflag == 2) {
      cout << "condition # estimate      relative residual" 
	   << "   iteration" << endl;
      cout << setiosflags(ios::scientific | ios::uppercase);
      for (i=0; i<num_iter; i++) {
	cout << " " 
	     << setw(17) << setprecision(10) << 0
	     << "       " 
	     << setw(17) << setprecision(10) << fabs(norms[i])
	     << "        " 
	     << i+1 << endl;
      }
    }
    cout << resetiosflags(ios::scientific);
    cout << resetiosflags(ios::uppercase);
    cout << setprecision(6);
  }
}

void CLOP_solver::apply_preconditioner(Epetra_Vector* r, Epetra_Vector* z)
{
  int i;
  double dprod;
  myzero(z_overlap, ndof_overlap);
  //
  // subdomain preconditioning with symmetric multiplicative correction
  // based on stored search directions
  //
  if ((pre_type_orthog == 1) && (n_orthog_used > 0)) {
    if (pre_type_coarse == 1) {
      /*
      PhiT->Multiply(false, *r, *PhiTr);
      PhiTr_gathered->Import(*PhiTr, *Importer_coarse, Insert);
      cout << *PhiTr_gathered << endl;
      */
      search_correction(r, z, 0.0);
      mat_vec_prod(z, ApSt_red);
      ApSt_red->Update(1.0, *r, -1.0);
      rOverlap->Import(*ApSt_red, *ImporterST2O, Insert);
      for (i=0; i<npart; i++) Asub[i].subpre(r_overlap, z_overlap, rhs_work, 
					     sol_work, tmp_work);
      vSt_red->Export(*zOverlap, *ExporterO2ST, Add);
      mat_vec_prod(vSt_red, gSt_red);
      ApSt_red->Update(-1.0, *gSt_red, 1.0);
      search_correction(ApSt_red, z, 1.0);
      z->Update(1.0, *vSt_red, 1.0);
      mat_vec_prod(z, ApSt_red);
      ApSt_red->Update(1.0, *r, -1.0);
      coarse_correction(ApSt_red, vSt_red);
      z->Update(1.0, *vSt_red, 1.0);
    }
    else {
      rOverlap->Import(*r, *ImporterST2O, Insert);
      for (i=0; i<npart; i++) Asub[i].subpre(r_overlap, z_overlap, rhs_work, 
					     sol_work, tmp_work);
      z->Export(*zOverlap, *ExporterO2ST, Add);
      coarse_correction(r, vSt_red);
      z->Update(1.0, *vSt_red, 1.0);  
      mat_vec_prod(z, ApSt_red);
      ApSt_red->Update(1.0, *r, -1.0);
      search_correction(ApSt_red, z, 1.0);
    }
  }
  else {
    rOverlap->Import(*r, *ImporterST2O, Insert);
    for (i=0; i<npart; i++) Asub[i].subpre(r_overlap, z_overlap, rhs_work, 
					   sol_work, tmp_work);
    z->Export(*zOverlap, *ExporterO2ST, Add);
    if (pre_type_coarse == 1) {
      mat_vec_prod(z, ApSt_red);
      ApSt_red->Update(1.0, *r, -1.0);
      coarse_correction(ApSt_red, vSt_red);
    }
    else coarse_correction(r, vSt_red);
    z->Update(1.0, *vSt_red, 1.0);
  }
}

void CLOP_solver::two_steps_CGS(int gmres_iter, Epetra_Vector* r)
{
  int i;
  double dprod, normr;
  double *A = VV;
  double *X = r->Values();
  double *Y = gmres_vec;
  double *Z = gmres_sum;
  int M = r->MyLength();
  int N = gmres_iter+1;
  int LDA = M;
  char TRANS = 'T'; double ALPHA(1); double BETA(0);
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
  else myzero(gmres_vec, N);
  Comm.SumAll(gmres_vec, gmres_sum, N);
  for (i=0; i<N; i++) HH[i]  = gmres_sum[i];
  TRANS = 'N'; ALPHA = -1; BETA = 1;
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, Z, BETA, X);
  TRANS = 'T'; ALPHA =  1; BETA = 0;
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
  else myzero(gmres_vec, N);
  Comm.SumAll(gmres_vec, gmres_sum, N);
  for (i=0; i<N; i++) HH[i] += gmres_sum[i];
  TRANS = 'N'; ALPHA = -1; BETA = 1;
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, Z, BETA, X);
  r->Dot(*r, &dprod);
  normr = sqrt(dprod);
  HH[N] = normr;
  int ibeg = M*N;
  double *V = &VV[ibeg];
  for (i=0; i<M; i++) V[i] = X[i]/normr;
  for (i=0; i<=N; i++) zz[i] = -HH[i];
  /*  
  TRANS = 'T'; ALPHA = 1; BETA = 0;
  for (i=0; i<N+1; i++) {
    EB.GEMV(TRANS, M, N+1, ALPHA, A, LDA, &VV[M*i], BETA, Y);
    Comm.SumAll(gmres_vec, gmres_sum, N+1);
    if (MyPID == 0) {
      cout << "orthogonality for i = " << i << endl;
      for (int j=0; j<N+1; j++) cout << gmres_sum[j] << endl;
    }
  }
  */
}

void CLOP_solver::hessenberg_qr(int gmres_iter)
{
  int i, ibeg;
  double w1, w2, vsign, normz;
  for (i=0; i<gmres_iter; i++) {
    w1 = cc[i]*zz[i] - ss[i]*zz[i+1];
    w2 = ss[i]*zz[i] + cc[i]*zz[i+1];
    zz[i] = w1;
    zz[i+1] = w2;
  }
  ibeg = gmres_iter*maxiter;
  for (i=0; i<gmres_iter; i++) RR[ibeg+i] = -zz[i];
  i = gmres_iter;
  gmres_givens(zz[i], zz[i+1], cc[i], ss[i]);
  if (ss[i] > 0) vsign = 1;
  else vsign = -1;
  normz = sqrt(zz[i]*zz[i] + zz[i+1]*zz[i+1]);
  RR[ibeg+gmres_iter] = -normz*vsign;
  if (i == 0) norms[i] = ss[i];
  else norms[i] = norms[i-1]*ss[i];
}  

void CLOP_solver::gmres_givens(double a, double b, double & c, double & s)
{
  double tau;
  if (b == 0) {
    c = 1;
    s = 0;
  }
  else { 
    if (fabs(a) > fabs(b)) {
      tau = b/a;
      c = 1/sqrt(1 + tau*tau);
      s = -tau*c;
    }
    else {
      tau = a/b;
      s = 1/sqrt(1 + tau*tau);
      c = -tau*s;
    }
  }
}

void CLOP_solver::construct_solution(int gmres_iter, double normb)
{
  int i, j, ii, n;
  double w1, w2;
  myzero(zz, gmres_iter+1);
  zz[0] = 1;
  for (i=0; i<gmres_iter; i++) {
    w1 = cc[i]*zz[i] - ss[i]*zz[i+1];
    w2 = ss[i]*zz[i] + cc[i]*zz[i+1];
    zz[i] = w1;
    zz[i+1] = w2;
  }
  n = gmres_iter;
  for (i=0; i<n; i++) {
    ii = n - i - 1; 
    for (j=ii+1; j<n; j++) zz[ii] -= RR[maxiter*j+ii]*zz[j];
    zz[ii] /= RR[maxiter*ii+ii];
  }
  char TRANS = 'N';
  double *A = VV;
  double *X = zz;
  double *Y = rSt_red->Values();
  int M = rSt_red->MyLength();
  int N = gmres_iter;
  int LDA = M;
  double ALPHA(normb); double BETA(0);
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
  apply_preconditioner(rSt_red, uSt_red);
}

void CLOP_solver::search_correction(Epetra_Vector *r, Epetra_Vector *u, 
				    double BETA2)
{
  char TRANS = 'T';
  double *A = P_matrix->Values();
  double *X = r->Values();
  double *Y = ortho_vec;
  double *Z = u->Values();
  int M = P_matrix->MyLength();
  int N = n_orthog_used;
  int LDA = M;
  double ALPHA(1); double BETA(0);
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
  else myzero(Y, N);
  Comm.SumAll(ortho_vec, ortho_sum, N);
  if (orthog_option == 2) {
    TRANS = 'N'; int INFO;
    EL.GETRS(TRANS, N, 1, PAP_sum, max_orthog, IPIV, ortho_sum, N, &INFO);
  }
  memcpy(Y, ortho_sum, N*sizeof(double));
  //  for (int i=0; i<N; i++) Y[i] /= pAp_vec[i];
  TRANS = 'N';
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, Y, BETA2, Z);
}

void CLOP_solver::remove_projection_search(Epetra_Vector* v)
{
  int j, number_cgs(2);
  char TRANS;
  double ALPHA, BETA;
  double *A = AP_matrix->Values();
  double *P =  P_matrix->Values();
  double *X = v->Values();
  double *Y = ortho_vec;
  int M = AP_matrix->MyLength();
  int N = n_orthog_used;
  int LDA = M;
  if (orthog_option == 1) {
    for (j=0; j<number_cgs; j++) {
      TRANS = 'T'; ALPHA = 1; BETA = 0;
      if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
      else myzero(Y, N);
      Comm.SumAll(ortho_vec, ortho_sum, N);
      memcpy(Y, ortho_sum, N*sizeof(double));
      //    for (int i=0; i<N; i++) Y[i] = ortho_sum[i]/pAp_vec[i];
      TRANS = 'N'; ALPHA = -1; BETA = 1;
      if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, P, LDA, Y, BETA, X);
    }
    /*
      TRANS = 'T'; ALPHA = 1; BETA = 0;
      EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
      Comm.SumAll(ortho_vec, ortho_sum, N);
      if (MyPID == 0) {
      if (MyPID == 0) {
      cout << "orthog_sum_after = " << endl;
      for (int i=0; i<n_orthog_used; i++) cout << ortho_sum[i] << endl;
      }
      }
    */
  }
  if (orthog_option == 2) {
    TRANS = 'T'; ALPHA = 1; BETA = 0;
    if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
    else myzero(Y, N);
    Comm.SumAll(ortho_vec, ortho_sum, N);
    TRANS = 'N'; int INFO;
    EL.GETRS(TRANS, N, 1, PAP_sum, max_orthog, IPIV, ortho_sum, N, &INFO);
    memcpy(Y, ortho_sum, N*sizeof(double));
    TRANS = 'N'; ALPHA = -1; BETA = 1;
    if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, P, LDA, Y, BETA, X);
  }
}

void CLOP_solver::store_search(double pAp)
{
  int i, j, INFO;
  if (n_orthog_used < max_orthog) {
    double *P  =  P_matrix->Values();
    double *AP = AP_matrix->Values();
    double *pvec  =  pSt_red->Values();
    double *Apvec = ApSt_red->Values();
    int M = P_matrix->MyLength();
    int ibeg =  M*n_orthog_used;
    memcpy( &P[ibeg],  pvec, M*sizeof(double));
    memcpy(&AP[ibeg], Apvec, M*sizeof(double));
    //    pAp_vec[n_orthog_used] = pAp;
    pAp_vec[n_orthog_used] = 1;
    double ALPHA = 1/sqrt(pAp);
    EB.SCAL(M, ALPHA, &P[ibeg]);
    EB.SCAL(M, ALPHA, &AP[ibeg]);
    n_orthog_used++;
    if (orthog_option == 2) {
      char TRANSA('T'), TRANSB('N');
      ALPHA = 1; double BETA(0);
      int N = n_orthog_used; int Nt2 = N*2; int Nt3 = N*3;
      int ibeg = M*(N - 1);
      if (M > 0) {
	EB.GEMV(TRANSA, M, N, ALPHA,        P, M, &AP[ibeg], BETA, PAP);
	EB.GEMV(TRANSA, M, N, ALPHA,       AP, M,  &P[ibeg], BETA, &PAP[N]);  
      }
      else myzero(PAP, Nt2);
      Comm.SumAll(PAP, &PAP[Nt2], Nt2);
      ibeg = max_orthog*(N - 1);
      for (i=0; i<n_orthog_used; i++) {
	PAP_store[(N-1)+i*max_orthog] = PAP[Nt3+i];
	PAP_store[i+ibeg]             = PAP[Nt2+i];
      }
      for (i=0; i<n_orthog_used; i++) {
	for (j=0; j<n_orthog_used; j++) {
	  int ii = i+j*max_orthog;
	  PAP_sum[ii] = PAP_store[ii];
	}
      }
      EL.GETRF(N, N, PAP_sum, max_orthog, IPIV, &INFO);
    }
  }
}

void CLOP_solver::coarse_correction(const Epetra_Vector* r, Epetra_Vector* u)
{
  if (ncdof > ndof_global_red) {
    u->PutScalar(0.0);
    return;
  }
  int nrhs(1);
  //
  // calculate Phi^T * r
  //
  PhiT->Multiply(false, *r, *PhiTr);
  //
  // gather PhiTr into PhiTr_gathered
  //
  PhiTr_gathered->Import(*PhiTr, *Importer_coarse, Insert);
  //
  // solve coarse problem
  //
  if (MyPID == 0) Kc_fac->sol(nrhs, rhs_coarse, sol_coarse, temp_coarse);
  //
  // scatter coarse problem solution
  //
  PhiTr->Export(*CSol_gathered, *Importer_coarse, Insert);
  //
  // calculate Phi * CSol
  //
  Phi->Multiply(false, *PhiTr, *u);
}

void CLOP_solver::calculate_multipliers(Epetra_Vector* uStand, 
                double & norm_rconstraint, double & norm_conerror)
{
  int i, j, NumEntries, *Indices;
  double *Values, rconstraint;
  vStand->ExtractView(&Values);
  for (i=0; i<nmycon; i++) {
    lambda_local[i] = 0;
    if (mycdof[i] >= 0) lambda_local[i] = Values[mycdof[i]];
  }
  /*
  char fname[101];
  sprintf(fname,"%s%d","test", MyPID);
  sprintf(fname, "%s.dat", fname);
  ofstream ffout;
  ffout.open(fname);
  int N;
  N = Exporter_lam->SourceMap().NumMyElements();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << Exporter_lam->SourceMap().GID(i) << endl;

  N = Exporter_lam->TargetMap().NumMyElements();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << Exporter_lam->TargetMap().GID(i) << endl;

  N = ConStandard->NumMyRows();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << ConStandard->GRID(i) << endl;

  N = ConStandard->NumMyCols();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << ConStandard->GCID(i) << endl;

  N = CtT->NumMyCols();
  ffout << N << endl;
  for (i=0; i<N; i++) ffout << CtT->GCID(i) << endl;
  
  ffout << ConStandard->MaxNumEntries() << endl;
  for (i=0; i<ConStandard->NumMyRows(); i++) {
    ConStandard->ExtractMyRowView(i, NumEntries, Values, Indices);
    ffout << NumEntries << endl;
    for (j=0; j<NumEntries; j++) ffout << Indices[j] << " " 
				      << Values[j] << endl;
  }

  ffout << CtT->MaxNumEntries() << endl;
  for (i=0; i<CtT->NumMyRows(); i++) {
    CtT->ExtractMyRowView(i, NumEntries, Values, Indices);
    ffout << NumEntries << endl;
    for (j=0; j<NumEntries; j++) ffout << Indices[j] << " " 
				      << Values[j] << endl;
  }

  ffout.close();
  */
  Lambda->Export(*Lambda_local, *Exporter_lam, Insert);
  CtT->Multiply(false, *Lambda, *wStand);
  //  cout << *wStand << endl;
  vStand->Update(-1.0, *wStand, 1.0); 
  vStand->Dot(*vStand, &norm_rconstraint);
  norm_rconstraint = sqrt(norm_rconstraint);
  ConStandard->Multiply(true, *uStand, *ConError);
  double norm_uStand, norm_ConError, norm_ConStandard;
  uStand->NormInf(&norm_uStand);
  ConError->NormInf(&norm_ConError);
  norm_ConStandard = ConStandard->NormInf();
  norm_conerror = norm_ConError/norm_ConStandard/norm_uStand;
  //  Lambda->Export(*Lambda_local, *Exporter_lam, Insert);
}

void CLOP_solver::mpcforces(Epetra_Vector* uLocal, 
			    Epetra_Import* ImporterStLam2Loc)
{
  if (ncon_global == 0) return;
  uLocal->Import(*Lambda, *ImporterStLam2Loc, Insert);
}

void CLOP_solver::mat_vec_prod(Epetra_Vector* u, Epetra_Vector* Au)
{
  ASt_red->Multiply(false, *u, *Au);
}

void CLOP_solver::calculate_condition(int miter)
{
  int i, j, INFO, ip1, one = 1;
  double Z, WORK;
  char N = 'N';
  if (miter == 1) {
    econa[0] = 1;
    return;
  }
  for (i=0; i<miter; i++) {
    ip1 = i + 1;
    Dtri[0] = pApa[0]/rhoa[0]/rhoa[0];
    for (j=1; j<ip1; j++) {
      Dtri[j]   = (pApa[j-1]*betaa[j]*betaa[j]+pApa[j])/rhoa[j]/rhoa[j];
      Etri[j-1] = -pApa[j-1]*betaa[j]/rhoa[j-1]/rhoa[j];
    }
    dstev_(&N, &ip1, Dtri, Etri, &Z, &one, &WORK, &INFO, 1); 
    if (INFO != 0) {
      cout << "error in call to DSTEV in CLASP_solver::calculate_condition" 
	   << endl;
      cout << "INFO = " << INFO << endl;
    }
    econa[i] = Dtri[i]/Dtri[0];
  }
}

void CLOP_solver::spmat_datfile(const Epetra_CrsMatrix & A, char fname[], 
				int opt)
{
  int i, j, NumEntries, *Indices, grow, gcol;
  double *Values;
  ofstream ffout;
  ffout.open(fname);
  for (i=0; i<A.NumMyRows(); i++) { 
    A.ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (opt == 1)
	ffout << i+1 << " " << Indices[j]+1 << setw(22) << setprecision(15)
	     << Values[j] << endl;
      if (opt == 2) {
	grow = A.GRID(i); gcol = A.GCID(Indices[j]);
	ffout << grow+1 << " " << gcol+1 << setw(22) 
	     << setprecision(15) << Values[j] << endl;
      }
    }
  }
  ffout.close();
}
