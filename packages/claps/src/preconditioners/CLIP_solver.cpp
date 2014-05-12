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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <assert.h>
//#include <climits>
#include "CLIP_solver.hpp"
#include "myzero.hpp"
#include <algorithm>

#include "Claps_ConfigDefs.hpp"  // for definition of F77_FUNC

#define DSTEV_F77 F77_FUNC(dstev,DSTEV)
//#define DSTEV_F77 dstev_

extern "C"{
  void metis_partgraphrecursive(int *n, int *xadj, int *adjncy, int *vwgt,
				int *adjwgt, int *wgtflag, int *numflag, 
				int *nparts, int *options, int *edgecut, 
				int *part);
  void metis_partgraphkway(int *n, int *xadj, int *adjncy, int *vwgt,
			   int *adjwgt, int *wgtflag, int *numflag, 
			   int *nparts, int *options, int *edgecut, int *part);
  void DSTEV_F77(char* N, int* i, double Dtri[], double Etri[], double* Z, 
		 int* one, double* WORK, int* INFO, long lengthN); 
}

CLIP_solver::CLIP_solver(const Epetra_CrsMatrix* ASub_,
			 const Epetra_IntVector* NodalDofs_,
			 const Epetra_MultiVector* Coords_,
			 const Epetra_CrsMatrix* ConStandard_,
			 const double* clip_params_)
  : ASub(ASub_), NodalDofs(NodalDofs_), Coords(Coords_), 
    ConStandard(ConStandard_), clip_params(clip_params_), Comm(ASub->Comm())
{

  zero_pointers();

  cdof_option         = int(clip_params[0]); 
  solver_tol          = clip_params[1];
  maxiter             = int(clip_params[2]);
  max_orthog          = int(clip_params[3]);
  atype               = int(clip_params[4]);
  ndim                = int(clip_params[5]);
  local_solver        = int(clip_params[6]);
  prt_debug           = int(clip_params[7]);
  prt_summary         = int(clip_params[8]);
  chk_sub_singularity = int(clip_params[9]);
  krylov_method       = int(clip_params[10]);
  scale_option        = int(clip_params[11]);
  num_rigid_mode      = int(clip_params[12]);

  assert (Coords->NumVectors() == 3);
  assert (Coords->ConstantStride() == true);

  double starttime, endtime;
  Comm.Barrier();
  MyPID = Comm.MyPID();
  if (MyPID == 0) starttime = MPI_Wtime();
  NumProc = Comm.NumProc();
  ndof_global = ConStandard->NumGlobalRows();
  print_flag = -1;
  if (MyPID == 0) {
    print_flag = prt_summary + 10*prt_debug;
    fout.open("CLIP_solver.data");
    fout << "----------------- CLIP solver summary information "
	 << "-----------------" << std::endl;
  }
  if (print_flag > 0) fout << "number of global dofs        = " << ndof_global 
			   << std::endl;
  //
  // process constraint equations
  //
  process_constraints();
  Exporter = new Epetra_Export(*SubMap, *StandardMap);
  Importer = new Epetra_Import(*SubMap, *StandardMap);
  ndof_sub = SubMap->NumMyElements();
  //
  // determine components for subdomain
  //
  determine_components();
  //
  // determine dof sets
  //
  determine_dof_sets();
  //
  // determine corner dofs
  //
  determine_corner_dofs();
  //
  // determine if extra corners are needed to remove singularities
  //
  if (chk_sub_singularity == 1) determine_extra_corners();
  //
  // modify dof sets to account for selected corners and which types
  // of constraints are selected
  //
  modify_dof_sets();
  //
  // factor substructure matrices
  //
  factor_sub_matrices();
  //
  // calculate coarse interpolation matrices and coarse stiffness matrix
  //
  calculate_coarse();
  //
  // initialize for pcg iterations
  //
  pcg_init();
  //
  //  const Epetra_MpiComm &empicomm = 
  //               dynamic_cast<const Epetra_MpiComm &>(Comm);
  //  mpicomm = empicomm.Comm();
  Comm.Barrier();
  if (print_flag > 0) {
    endtime = MPI_Wtime();
    fout << "elapsed time for clip solver init  = " 
	 << endtime-starttime << " seconds" << std::endl;
    if (krylov_method != 1) fout << "pcg solver will be used" << std::endl;
    if (krylov_method == 1) fout << "gmres solver will be used" << std::endl;
    fout.close();
  }
}

CLIP_solver::~CLIP_solver()
{
  delete Importer; delete ImporterB; delete ImporterSt2B; 
  delete Importer_Kc; delete Exporter; delete ExporterB; delete Exporter_lam; 
  delete Exporter_Kc; delete ImporterStand; delete ExporterStand; 
  delete Lambda; delete Lambda_local; delete ConError; delete ApB_Sub; 
  delete ApB_St; delete pB_Sub; delete pB_St; delete zB_Sub; delete zB_St; 
  delete uB_St; delete u_St; delete vStand; delete r_St; delete r_Sub; 
  delete rB_St; delete rB_Sub; delete work_Sub; delete work_Kc; delete rBSub; 
  delete workB_St; delete work_St; delete work_St2; delete vec_PhiB; 
  delete prod_PhiB; delete vec1_ASub; delete vec2_ASub; delete wStand; 
  delete StandardMap; delete SubMap; delete RowMapMyCon; delete MapB_Sub; 
  delete MapB_St; delete CtT; delete Tran; delete Block_Stiff; delete PhiB; 
  delete Coords_red; delete NodalDofs_red; delete AR; delete AI; 
  delete AKc; delete [] comp1; delete [] comp2; delete [] sub1; delete [] sub2;
  delete [] dset1; delete [] dset2; delete [] corner_flag; delete [] mycdof; 
  delete [] dofI; delete [] dofB; delete [] dofC; delete [] dofR; 
  delete [] lambda; delete [] lambda_local; delete [] weight; 
  delete [] ARinvCT; delete [] CARinvCT; delete [] lambda_e; delete [] RHS_cg; 
  delete [] SOL_cg; delete [] TEMP_cg; delete [] SOL_Kc; delete [] TEMP_Kc; 
  delete [] rcurra; delete [] rhoa; delete [] betaa; delete [] pApa; 
  delete [] Dtri; delete [] Etri; delete [] econa; delete [] P_ortho; 
  delete [] AP_ortho; delete [] orth1; delete [] orth2; delete [] owner_flag;
  delete [] VV; delete [] RR; delete [] HH; delete [] zz; delete [] cc;
  delete [] ss; delete [] norms; delete [] gmres_vec; delete [] gmres_sum;
  delete [] tied_down; delete Phir_St; delete [] PhirTPhir; delete Rhs_null;
  zero_pointers();
}

void CLIP_solver::process_constraints()
{
  int i, ScanSums, ndof_sub_all, *gdof_block;
  int NumEntries, *Indices;
  double *Values;
  //
  // construct block stiffness matrix Block_Stiff
  //
  ncon_global = 0;
  ndof_sub = ASub->NumMyRows();
  Comm.ScanSum(&ndof_sub, &ScanSums, 1);
  Comm.SumAll(&ndof_sub, &ndof_sub_all, 1);
  gdof_block = new int[ndof_sub];
  for (i=0; i<ndof_sub; i++) gdof_block[i] = ScanSums - ndof_sub + i;
  Epetra_Map Row_Map(ndof_sub_all, ndof_sub, gdof_block, 0, Comm);
  delete [] gdof_block;
  int *NumEntriesPerRow = new int[ndof_sub];
  for (i=0; i<ndof_sub; i++) {
    ASub->ExtractMyRowView(i, NumEntries, Values, Indices);
    NumEntriesPerRow[i] = NumEntries;
  }
  Block_Stiff = new Epetra_CrsMatrix(Copy, Row_Map, Row_Map,
				     NumEntriesPerRow, true);
  for (i=0; i<ndof_sub; i++) {
    ASub->ExtractMyRowView(i, NumEntries, Values, Indices);
    Block_Stiff->InsertMyValues(i, NumEntries, Values, Indices);
  }
  Block_Stiff->FillComplete();
  Block_Stiff->OptimizeStorage();
  double nnz_before = Block_Stiff->NumGlobalNonzeros();
  //
  // process constraints
  //
  ncon_global = ConStandard->NumGlobalCols();
  if (print_flag > 0) {
    fout << "number of global constraints = " << ncon_global << std::endl;
    fout << "constraint data ------------------------------------" << std::endl;
  }
 if (ncon_global > 0) {
    int nsub_gdofs, *sub_gdofs, flag(0), nx2, nx2_global, *x2_dof;
    nsub_gdofs = ASub->RowMap().NumMyPoints();
    sub_gdofs = new int[nsub_gdofs];
    ASub->RowMap().MyGlobalElements(sub_gdofs);
    //
    // transform constraints to standard form
    //
    CLOP_constraint *A;
    A = new CLOP_constraint(ConStandard, &fout, print_flag);
    flag = A->factor();
    A->Tran(Tran, RowMapMyCon, mycdof, nx2, x2_dof, nsub_gdofs, 
	    sub_gdofs, CtT);
    delete [] sub_gdofs;
    nmycon = RowMapMyCon->NumMyElements();
    int nn = ConStandard->DomainMap().NumMyPoints();
    lambda_local = new double[nmycon];
    lambda  = new double[nn];
    Lambda_local = new Epetra_Vector(View, *RowMapMyCon, lambda_local);
    Lambda       = new Epetra_Vector(View, ConStandard->DomainMap(), lambda);
    ConError     = new Epetra_Vector(ConStandard->DomainMap());
    Exporter_lam = new Epetra_Export(*RowMapMyCon, ConStandard->DomainMap());
    Comm.SumAll(&nx2, &nx2_global, 1);
    delete A;

    assert(nx2_global == 0);
    delete [] x2_dof;

    Epetra_CrsMatrix Tran_Sub(Copy, ASub->RowMap(), 0);
    Epetra_Import Import_Sub(ASub->RowMap(), Tran->RangeMap());
    Tran_Sub.Import(*Tran, Import_Sub, Insert);
    Tran_Sub.FillComplete(Tran->DomainMap(), Tran->RangeMap());

    copy_map(Tran->DomainMap(), StandardMap);
    copy_map(Tran_Sub.ColMap(), SubMap);
    /*
    char fname[101];
    sprintf(fname,"%s%d","test", MyPID);
    sprintf(fname, "%s.dat", fname);
    std::ofstream ffout;
    ffout.open(fname);
    ffout << ASub->NumMyRows() << std::endl;
    for (i=0; i<ASub->NumMyRows(); i++) ffout << ASub->GRID(i) << std::endl;
    ffout << StandardMap->NumMyElements() << std::endl;
    for (i=0; i<StandardMap->NumMyElements(); i++) 
      ffout << StandardMap->GID(i) << std::endl;
    ffout.close();
    */
    // replaced *StandardMap with ConStandard->RowMap() 
    Epetra_Export Exporter_temp(ASub->RowMap(), ConStandard->RowMap());
    Epetra_MultiVector  Coords_temp(ConStandard->RowMap(), 
				    Coords->NumVectors());
    Epetra_IntVector NodalDofs_temp(ConStandard->RowMap());
    Coords_temp.Export(*Coords, Exporter_temp, Insert);
    NodalDofs_temp.Export(*NodalDofs, Exporter_temp, Insert);
    Epetra_Import Importer_temp(*SubMap, ConStandard->RowMap());
    Coords_red = new Epetra_MultiVector( *SubMap, Coords->NumVectors());
    NodalDofs_red = new Epetra_IntVector(*SubMap);
    Coords_red->Import(Coords_temp, Importer_temp, Insert);
    NodalDofs_red->Import(NodalDofs_temp, Importer_temp, Insert);

    int ncol_Tran, ncol_Tran_all;
    ncol_Tran = Tran_Sub.NumMyCols();
    Comm.ScanSum(&ncol_Tran, &ScanSums, 1);
    Comm.SumAll(&ncol_Tran, &ncol_Tran_all, 1);
    gdof_block = new int[ncol_Tran];
    for (i=0; i<ncol_Tran; i++) gdof_block[i] = ScanSums - ncol_Tran + i;
    Epetra_Map Col_Map(ncol_Tran_all, ncol_Tran, gdof_block, 0, Comm);
    delete [] gdof_block;
    for (i=0; i<ndof_sub; i++) {
      Tran_Sub.ExtractMyRowView(i, NumEntries, Values, Indices);
      NumEntriesPerRow[i] = NumEntries;
    }
    Epetra_CrsMatrix *Block_Tran, *Block_TranT, *Block_KT;
    Block_Tran = new Epetra_CrsMatrix (Copy, Row_Map, Col_Map,
				       NumEntriesPerRow, true);
    for (i=0; i<ndof_sub; i++) {
      Tran_Sub.ExtractMyRowView(i, NumEntries, Values, Indices);
      Block_Tran->InsertMyValues(i, NumEntries, Values, Indices);
    }
    Block_Tran->FillComplete(Col_Map, Row_Map);
    Block_Tran->OptimizeStorage();
    construct_transpose(Block_Tran, Block_TranT);
    Block_TranT->OptimizeStorage();
    EpetraExtCD::MatrixMatrix::Multiply(*Block_Stiff, *Block_Tran, Block_KT);
    delete Block_Stiff;
    EpetraExtCD::MatrixMatrix::Multiply(*Block_TranT, *Block_KT, 
					Block_Stiff);
    Block_Stiff->OptimizeStorage();
    delete Block_Tran; delete Block_TranT; delete Block_KT;
    /*
    Epetra_Vector v1(Col_Map);
    Epetra_Vector v2(Row_Map);
    v1.PutScalar(1.0);
    Block_Stiff->Multiply(false, v1, v2);
    std::cout << v2 << std::endl;
    */
  }
  else {
    copy_map(ConStandard->RowMap(), StandardMap);
    copy_map(ASub->RowMap(), SubMap);
    copy_multivector(Coords, ASub->RowMap(), Coords_red);
    copy_intvector(NodalDofs, ASub->RowMap(), NodalDofs_red);
  }
 delete [] NumEntriesPerRow;
  u_St = new Epetra_Vector(*StandardMap);
  if (ncon_global > 0) {
    ImporterStand = new Epetra_Import(ASub->RowMap(), ConStandard->RowMap());
    ExporterStand = new Epetra_Export(ASub->RowMap(), ConStandard->RowMap());
    vec1_ASub = new Epetra_Vector(ASub->RowMap());
    vec2_ASub = new Epetra_Vector(ASub->RowMap());
    vStand = new Epetra_Vector(ConStandard->RowMap());
    wStand = new Epetra_Vector(ConStandard->RowMap());
    u_St->Random();
    Tran->Multiply(false, *u_St, *vStand);
    ConStandard->Multiply(true, *vStand, *Lambda);
    double inf_norm_error; Lambda->NormInf(&inf_norm_error);
    double con_error_norm = inf_norm_error/ConStandard->NormInf();
    int max_nnz_row = Tran->GlobalMaxNumEntries();
    double inf_norm_Tran = Tran->NormInf();
    double nnz_after  =   Block_Stiff->NumGlobalNonzeros();
    if (MyPID == 0) {
      if (print_flag > 0) {
	fout << "ratio of nnzs after static condensation = "
	     << nnz_after/nnz_before << std::endl;
	fout << "normalized constraint error check       = " << con_error_norm 
	     << std::endl;
	fout << "maximum nnz in any row of T matrix      = " << max_nnz_row 
	     << std::endl;
	fout << "infinity norm of T matrix               = " << inf_norm_Tran 
	     << std::endl;
	fout << "----------------------------------------------------" << std::endl;
      }
    }
    assert(con_error_norm < 1e-10);
  }
}

void CLIP_solver::copy_map(const Epetra_Map A, Epetra_Map* & B)
{
  int i, NumMyElements, NumGlobalElements, *MyGlobalElements;
  NumMyElements = A.NumMyElements();
  NumGlobalElements = A.NumGlobalElements();
  MyGlobalElements = A.MyGlobalElements();
  B = new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,
		     0, Comm);
  return;
}

void CLIP_solver::copy_multivector(const Epetra_MultiVector* A,
				   const Epetra_Map RowMap,
				   Epetra_MultiVector* & B)
{
  int NumVectors, MyLDA;
  double *valA;
  NumVectors = A->NumVectors();
  A->ExtractView(&valA, &MyLDA);
  B = new Epetra_MultiVector(Copy, RowMap, valA, MyLDA, NumVectors);
}

void CLIP_solver::copy_intvector(const Epetra_IntVector* A,
				 const Epetra_Map RowMap,
				 Epetra_IntVector* & B)
{
  int *valA;
  A->ExtractView(&valA);
  B = new Epetra_IntVector(Copy, RowMap, valA);
}

void CLIP_solver::construct_transpose(Epetra_CrsMatrix* & A, 
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
  Epetra_Export Exporter_temp(B.RowMap(), A->DomainMap());
  AT->Export(B, Exporter_temp, Add);
  AT->FillComplete(A->RangeMap(), A->DomainMap());
  delete [] count; delete [] rowbeg; delete [] colidx; delete [] Bval;
}

void CLIP_solver::determine_components()
{
  int i, j, NumEntries, *Indices;
  double *Values;
  //
  // first determine dof connectivity matrix:
  //  a1[a2[i]:a2[i+1]-1] = dofs connected to dof i
  //
  int *a2 = new int[ndof_sub+1]; a2[0] = 0;
  int *arow, *iflag, narow, dof, dof2;
  arow = new int[ndof_sub];
  iflag = new int[ndof_sub]; myzero(iflag, ndof_sub);
  for (i=0; i<ndof_sub; i++) {
    narow = 0;
    Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      dof2 = Indices[j];
      if ((Values[j] != 0) && (iflag[dof2] == 0)) {
	arow[narow] = dof2;
	iflag[dof2] = 1;
	narow++;
      }
    }
    a2[i+1] = a2[i] + narow;
    for (j=0; j<narow; j++) iflag[arow[j]] = 0;
  }
  int *a1 = new int[a2[ndof_sub]];
  for (i=0; i<ndof_sub; i++) {
    narow = 0;
    Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      dof2 = Indices[j];
      if ((Values[j] != 0) && (iflag[dof2] == 0)) {
	arow[narow] = dof2;
	a1[a2[i] + narow] = dof2;
	iflag[dof2] = 1;
	narow++;
      }
    }
    for (j=0; j<narow; j++) iflag[arow[j]] = 0;
  }
  //
  // now determine subdomain components
  //  ncomp = number of components
  //  comp1[comp2[i]:comp2[i+1]-1] = subdomain dofs in component i
  //
  determine_components(a1, a2, ndof_sub, comp1, comp2, ncomp);
  //  std::cout << "MyPID, ncomp = " << MyPID << " " << ncomp << std::endl;
  delete [] a1; delete [] a2; delete [] arow; delete [] iflag;
}

void CLIP_solver::determine_dof_sets()
{
  int i, j, NumEntries, *Indices, dof, *loc_dofs;
  int ScanSums;
  NodalDofs->ExtractView(&loc_dofs);
  //
  // each subdomain component is now treated and numbered as an
  // individual subdomain
  //  ncomp_sum = total number of subdomains
  //
  Comm.ScanSum(&ncomp, &ScanSums, 1);
  Comm.SumAll(&ncomp, &ncomp_sum, 1);
  //  std::cout << "MyPID, ScanSums = " << MyPID << " " << ScanSums << std::endl;
  //
  // first determine subdomains containing each dof
  //  sub1[sub2[i]:sub2[i+1]-1] = subdomains containing dof i
  //
  int *columns = new int[ncomp];
  sub_begin = ScanSums - ncomp;
  for (i=0; i<ncomp; i++) columns[i] = sub_begin + i;
  Epetra_Map ColMap(ncomp_sum, ncomp, columns, 0, Comm);
  delete [] columns;
  Epetra_CrsGraph Subdomains(Copy, *SubMap, ColMap, 1, true);
  for (i=0; i<ncomp; i++) {
    for (j=comp2[i]; j<comp2[i+1]; j++) {
      dof = comp1[j];
      Subdomains.InsertMyIndices(dof, 1, &i);
    }
  }
  Subdomains.FillComplete(ColMap, *StandardMap);
  Epetra_CrsGraph SubStand(Copy, *StandardMap, 0);
  SubStand.Export(Subdomains, *Exporter, Insert);
  SubStand.FillComplete(ColMap, *StandardMap);
  Epetra_CrsGraph SubsForDofs(Copy, *SubMap, 0);
  SubsForDofs.Import(SubStand, *Importer, Insert);
  SubsForDofs.FillComplete(ColMap, *StandardMap);
  sub1 = new int[SubsForDofs.NumMyNonzeros()];
  sub2 = new int[ndof_sub+1]; sub2[0] = 0;
  for (i=0; i<ndof_sub; i++) {
    SubsForDofs.ExtractMyRowView(i, NumEntries, Indices);
    for (j=0; j<NumEntries; j++) {
      sub1[sub2[i]+j] = SubsForDofs.GCID(Indices[j]);
      sub2[i+1] = sub2[i] + NumEntries;
    }
    std::sort(&sub1[sub2[i]], &sub1[sub2[i+1]]);
  }
  /*
  if (MyPID == 3) {
    for (i=0; i<ndof_sub; i++) {
      for (j=sub2[i]; j<sub2[i+1]; j++) std::cout << sub1[j] << " ";
      std::cout << std::endl;
    }
  }
  */
  //
  // map subdomains for each dof to a double (subdomain numbers in each 
  // segment of sub1 assumed sorted)
  //
  int max_num_sub(0), num_sub, subnum;
  for (i=0; i<ndof_sub; i++) {
    num_sub = sub2[i+1] - sub2[i];
    if (num_sub > max_num_sub) max_num_sub = num_sub;
  }
  double *randa = new double[max_num_sub];
  for (i=0; i<max_num_sub; i++) randa[i] = 0.777*rand()/RAND_MAX;
  double *dval = new double[ndof_sub];
  double *dval_orig = new double[ndof_sub];
  double dsum;
  for (i=0; i<ndof_sub; i++) {
    dsum = 0;
    for (j=sub2[i]; j<sub2[i+1]; j++) {
      subnum = sub1[j];
      dsum += static_cast<double>(subnum + 1) / ncomp_sum * randa[j-sub2[i]] +
	loc_dofs[i] * 0.7777777;
    }
    dval[i] = dsum;
    dval_orig[i] = dsum;
  }
  delete [] randa;
  //
  // next sort mapped values and determine number of dof sets
  //
  std::sort(&dval[0], &dval[ndof_sub]);
  ndof_set = 0; double prev_val = -7;
  for (i=0; i<ndof_sub; i++) {
    if (dval[i] != prev_val) {
      prev_val = dval[i];
      ndof_set++;
    }
  }
  double *values = new double[ndof_set];
  ndof_set = 0; prev_val = -7;
  for (i=0; i<ndof_sub; i++) {
    if (dval[i] != prev_val) {
      prev_val = dval[i];
      values[ndof_set] = dval[i];
      ndof_set++;
    }
  }
  //
  // next determine dof set number for each subdomain dof
  //  ival[i] = dof set for subdomain dof i
  //
  int *ival = new int[ndof_sub];
  int *icount = new int[ndof_set]; myzero(icount, ndof_set);
  for (i=0; i<ndof_sub; i++) {
    double *lower = std::lower_bound(&values[0], &values[ndof_set], dval_orig[i]);
    ival[i] = (int) (lower - &values[0]);
    icount[ival[i]]++;
    /*
    if (MyPID == 0) std::cout << "val, ival = " << dval_orig[i] << " " 
    			 << ival[i] << std::endl;
    */
  }
  delete [] dval; delete [] dval_orig; delete [] values;
  //
  // determine nodes in each dof set
  //  dset1[dset2[i]:dset2[i+1]-1] = subdomain dofs in dof set i
  //
  dset1 = new int[ndof_sub];
  dset2 = new int[ndof_set+1]; dset2[0] = 0;
  for (i=0; i<ndof_set; i++) dset2[i+1] = dset2[i] + icount[i];
  myzero(icount, ndof_set);
  for (i=0; i<ndof_sub; i++) {
    dset1[dset2[ival[i]] + icount[ival[i]]] = i;
    icount[ival[i]]++;
  }
  delete [] ival; delete [] icount;
  //
  // next determine ownership for each subdomain dof
  //  owner_flag = true if subdomain owns subdomain dof i, false otherwise
  //
  owner_flag = new bool[ndof_sub];
  for (i=0; i<ndof_sub; i++) owner_flag[i] = false;
  int base = ncomp_sum + 7;
  for (i=0; i<ndof_sub; i++) {
    int gmin = base;
    for (j=sub2[i]; j<sub2[i+1]; j++) if (sub1[j] < gmin) gmin = sub1[j];
    if ((gmin >= sub_begin) && (gmin < (sub_begin + ncomp))) owner_flag[i] = true;
  }

  /*
  int *gdofs, MyLDA;
  gdofs = SubMap->MyGlobalElements();
  double *xyz;
  Coords_red->ExtractView(&xyz, &MyLDA);
  std::cout << *ConStandard << std::endl;
    if (MyPID == 1) {
    for (i=0; i<ndof_set; i++) {
      std::cout << "local, global, coords sub dofs in dof set: subs " << i << " :";
      dof = dset1[dset2[i]];
      for (j=sub2[dof]; j<sub2[dof+1]; j++) std::cout << sub1[j] << " ";
      std::cout << std::endl;
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	std::cout << dset1[j] << " " << gdofs[dset1[j]] << ": " << xyz[dset1[j]] 
	     << " " << xyz[dset1[j]+MyLDA] << " " << xyz[dset1[j]+2*MyLDA]
	     << std::endl;
      }
    }
  }
  */
  //
  // finally, for purposes of consistency between subdomains, order entries
  // in dset1 segments by increasing global dof value
  //
  /*
  int max_num_dof(0), num_dof;
  for (i=0; i<ndof_set; i++) {
    dof = dset1[dset2[i]];
    num_sub = sub2[dof+1] - sub2[dof];
    num_dof = dset2[i+1] - dset2[i];
    if ((num_sub > 1) && (num_dof > max_num_dof)) max_num_dof = num_dof;
  }
  int *dset_orig, *gdof_sort, sorted_loc, gdof, iflag;
  dset_orig = new int[max_num_dof];
  gdof_sort = new int[max_num_dof];
  for (i=0; i<ndof_set; i++) {
    dof = dset1[dset2[i]];
    num_sub = sub2[dof+1] - sub2[dof];
    if (num_sub > 1) {
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	dof = dset1[j];
	dset_orig[j-dset2[i]] = dof;
	gdof_sort[j-dset2[i]] = gdofs[dof];
      }
      num_dof = dset2[i+1] - dset2[i];
      std::sort(&gdof_sort[0], &gdof_sort[num_dof]);
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	dof = dset_orig[j-dset2[i]];
	gdof = gdofs[dof];
	sorted_loc = CRD_utils::find_index(&gdof_sort[0], num_dof, gdof);
	dset1[dset2[i]+sorted_loc] = dof;
      }
      iflag = 0;
      for (j=dset2[i]+1; j<dset2[i+1]; j++)
	if (gdofs[dset1[j-1]] >= gdofs[dset1[j]]) iflag = 1;
      assert (iflag == 0);
    }
  }
  delete [] dset_orig; delete [] gdof_sort;
  */
  /*
  if (MyPID == 0) {
    for (i=0; i<ndof_set; i++) {
      std::cout << "global sub dofs in dof set " << i << std::endl;
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	std::cout << gdofs[dset1[j]] << std::endl;
      }
    }
  }
  */  
}

void CLIP_solver::determine_corner_dofs()
{
  int i, j, k, num_sub, num_entries, dof, MyLDA, nexist, i1, i2, i3, *gdofs;
  double x, y, z;  
  gdofs = SubMap->MyGlobalElements();
  //
  // first pick singleton dof sets as corner dofs
  //
  corner_flag = new int[ndof_sub]; myzero(corner_flag, ndof_sub);
  for (i=0; i<ndof_set; i++) {
    num_entries = dset2[i+1] - dset2[i];
    dof = dset1[dset2[i]];
    if (num_entries == 1) corner_flag[dof] = 1;
  }
  //
  // next pick corner dofs for edges
  //
  int rot_flag;
  int *adof = new int[ndof_sub];
  bool *adof_flag = new bool[ndof_sub];
  for (i=0; i<ndof_sub; i++) adof_flag[i] = false;
  int ndof_cand;
  double *xyz;
  Coords_red->ExtractView(&xyz, &MyLDA);
  for (i=0; i<ndof_set; i++) {
    num_entries = dset2[i+1] - dset2[i];
    dof = dset1[dset2[i]];
    num_sub = sub2[dof+1] - sub2[dof];
    rot_flag = determine_shell(i);
    if ((num_entries > 1) && (num_sub > 1)) {
      if ((num_sub > 2) || (ndim == 2) || (rot_flag == 1)) {
	get_cand_dofs(i, adof, adof_flag, ndof_cand);
	nexist = 0;
	if (MyPID == -1) {
	  std::cout << "num_sub  = " << num_sub << std::endl;
	  std::cout << "ndim     = " << ndim << std::endl;
	  std::cout << "rot_flag = " << rot_flag << std::endl;
	  std::cout << "coordinates for edge dof set" << std::endl;
	  for (j=dset2[i]; j<dset2[i+1]; j++) {
	    dof = dset1[j];
	    std::cout << xyz[dof] << " " << xyz[dof+MyLDA] << " " 
		 << xyz[dof+2*MyLDA] << std::endl;
	  }
	  std::cout << "ndof_cand = " << ndof_cand << std::endl;
	  std::cout << "coordinates of candidate dofs" << std::endl;
	  for (j=0; j<ndof_cand; j++) {
	    dof = adof[j];
	    std::cout << xyz[dof] << " " << xyz[dof+MyLDA] << " " 
		 << xyz[dof+2*MyLDA] << std::endl;
	  }
	}
	for (j=0; j<ndof_cand; j++) {
	  dof = adof[j];
	  if (corner_flag[dof] == 1) {
	    i1 = dof;
	    nexist++;
	  }
	}
	if (nexist == 0) {
	  get_two_points(xyz, MyLDA, adof, ndof_cand, i1, i2);
	  corner_flag[i1] = 1;
	  corner_flag[i2] = 1;
	}
	if (nexist == 1) {
	  get_one_point( xyz, MyLDA, adof, ndof_cand, i1, i2);
	  corner_flag[i2] = 1;
	}
      }
    }
  }
  //
  // next pick corner dofs for faces
  //
  for (i=0; i<ndof_set; i++) {
    num_entries = dset2[i+1] - dset2[i];
    dof = dset1[dset2[i]];
    num_sub = sub2[dof+1] - sub2[dof];
    rot_flag = determine_shell(i);
    if ((num_entries > 1) && (num_sub > 1)) {
      if ((num_sub == 2) && (ndim == 3) && (rot_flag == 0)) {
	get_cand_dofs(i, adof, adof_flag, ndof_cand);
	if (MyPID == -1) {
	  std::cout << "coordinates for face dof set" << std::endl;
	  for (j=dset2[i]; j<dset2[i+1]; j++) {
	    dof = dset1[j];
	    std::cout << xyz[dof] << " " << xyz[dof+MyLDA] << " " 
		 << xyz[dof+2*MyLDA] << std::endl;
	  }
	  std::cout << "ndof_cand = " << ndof_cand << std::endl;
	  std::cout << "coordinates of candidate dofs" << std::endl;
	  for (j=0; j<ndof_cand; j++) {
	    dof = adof[j];
	    std::cout << xyz[dof] << " " << xyz[dof+MyLDA] << " " 
		 << xyz[dof+2*MyLDA] << std::endl;
	  }
	}
	nexist = 0;
	for (j=0; j<ndof_cand; j++) {
	  dof = adof[j];
	  if (corner_flag[dof] == 1) {
	    i1 = dof;
	    nexist++;
	  }
	}
	if (nexist == 0) {
	  get_first_point(adof, ndof_cand, i1);
	  corner_flag[i1] = 1;
	  get_one_point(xyz, MyLDA, adof, ndof_cand, i1, i2);
	  corner_flag[i2] = 1;
	  get_third_point(xyz, MyLDA, adof, ndof_cand, i1, i2, i3);
	  corner_flag[i3] = 1;
	}
	if (nexist == 1) {
	  get_one_point(xyz, MyLDA, adof, ndof_cand, i1, i2);
	  corner_flag[i2] = 1;
	  get_third_point(xyz, MyLDA, adof, ndof_cand, i1, i2, i3);
	  corner_flag[i3] = 1;
	}      
      }
    }
  }
  delete [] adof; delete [] adof_flag;
  //
  // share corner dofs between subdomains
  //
  share_corner_info();
  max_added_corner = 0;
  /*
  for (i=0; i<ndof_sub; i++)
    if (corner_flag[i] > 0) { 
      x = xyz[i]; y = xyz[i+MyLDA]; z = xyz[i+2*MyLDA];
      std::cout << "MyPID, corner, coords = " << MyPID << " " << gdofs[i] << " " 
	   << x << " " << y << " " << z << std::endl;
    }
  */
}

void CLIP_solver::share_corner_info()
{
  Epetra_IntVector C_sub(View, *SubMap, corner_flag);
  Epetra_IntVector C_all(*StandardMap);
  C_all.Export(C_sub, *Exporter, Add);
  C_sub.Import(C_all, *Importer, Insert);
  for (int i=0; i<ndof_sub; i++) if (corner_flag[i] > 0) corner_flag[i] = 1;
}

void CLIP_solver::determine_extra_corners()
{
  int i, *rowbeg, *colidx, *extra_corner, nextra, dof;
  double *vals, *Xvecs;
  //  myzero(corner_flag, ndof_sub);
  nR = 0; nC = 0;
  for (i=0; i<ndof_sub; i++) {
    if (corner_flag[i] == 0) nR++;
    else nC++;
  }
  dofR = new int[nR]; dofC = new int[nC];
  nR = 0; nC = 0;
  for (i=0; i<ndof_sub; i++) {
    if (corner_flag[i] == 0) {
      dofR[nR] = i;
      nR++;
    }
    else {
      dofC[nC] = i;
      nC++;
    }
  }
  bool *bound_flag = new bool[nR];
  for (i=0; i<nR; i++) {
    bound_flag[i] = false;
    dof = dofR[i];
    if ((sub2[dof+1] - sub2[dof]) > 1) bound_flag[i] = true;
  }
  gen_matrix(Block_Stiff, nR, dofR, rowbeg, colidx, vals);
  CRD_utils::subspace_iteration(nR, rowbeg, colidx, vals, bound_flag, 
				scale_option, nextra, extra_corner, AR, Xvecs);
  delete AR; delete [] dofR; delete [] dofC;
  for (i=0; i<nextra; i++) {
    dof = dofR[extra_corner[i]];
    corner_flag[dof] = 1;
    //    std::cout << "MyPID, nextra_corner = " << MyPID << " " << nextra << std::endl; 
  }
  Comm.MaxAll(&nextra, &max_added_corner, 1);
  delete [] bound_flag; delete [] rowbeg, delete [] colidx; delete [] vals;
  delete [] extra_corner;
  //
  // share corner dofs between subdomains
  //
  share_corner_info();
  /*
  double *xyz; int MyLDA;
  Coords_red->ExtractView(&xyz, &MyLDA);
  for (i=0; i<nextra; i++) {
    dof = adof[extra_corner[i]];
    std::cout << "MyPID, xyz extra = " << MyPID << " " << xyz[dof] << " " 
	 << xyz[dof+MyLDA] << " " << xyz[dof+2*MyLDA] << std::endl;
  }
  */
}

void CLIP_solver::modify_dof_sets()
{
  int i, j, dof, icount, nnz, nsub_for_dof, ndof_for_set;
  int ndof_set_orig, *dset2_orig;
  //
  // first remove existing corners from dof sets
  //
  dset2_orig = new int[ndof_set+1];
  for (i=0; i<=ndof_set; i++) dset2_orig[i] = dset2[i];
  ndof_set_orig = ndof_set;
  ndof_set = 0;
  nnz = 0;
  for (i=0; i<ndof_set_orig; i++) {
    icount = 0;
    for (j=dset2_orig[i]; j<dset2_orig[i+1]; j++) {
      dof = dset1[j];
      if (corner_flag[dof] == 0) {
	dset1[nnz] = dof;
	nnz++;
	icount++;
      }
    }
    if (icount > 0) {
      ndof_set++;
      dset2[ndof_set] = nnz;
    }
  }
  //
  // next remove dofs based on cdof_option
  //
  for (i=0; i<=ndof_set; i++) dset2_orig[i] = dset2[i];
  ndof_set_orig = ndof_set;
  ndof_set = 0;
  nnz = 0;
  for (i=0; i<ndof_set_orig; i++) {
    dof = dset1[dset2_orig[i]];
    nsub_for_dof = sub2[dof+1] - sub2[dof];
    ndof_for_set = dset2_orig[i+1] - dset2_orig[i];
    icount = 0;
    if ((cdof_option == 1) && (nsub_for_dof >  2)) icount = ndof_for_set;
    if ((cdof_option == 2) && (nsub_for_dof == 2)) icount = ndof_for_set;
    if ((cdof_option == 3) && (nsub_for_dof >  1)) icount = ndof_for_set;
    if (icount > 0) {
      for (j=dset2_orig[i]; j<dset2_orig[i+1]; j++) {
	dof = dset1[j];
	dset1[nnz] = dof;
	nnz++;
      }
      ndof_set++;
      dset2[ndof_set] = nnz;
    }
  }      
  delete [] dset2_orig;

  if (MyPID == -3) {
    std::cout << "dset2 = ";
    for (i=0; i<=ndof_set; i++) std::cout << dset2[i] << " ";
    std::cout << std::endl;
    double *xyz; int MyLDA;
    Coords_red->ExtractView(&xyz, &MyLDA);
    for (i=0; i<ndof_set; i++) {
      std::cout << "coordinates of dofs in dof_set " << i << std::endl;
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	dof = dset1[j];
	std::cout << xyz[dof] << " " << xyz[dof+MyLDA] << " " << xyz[dof+2*MyLDA]
	     << std::endl;
      }
    }
  }
}

void CLIP_solver::factor_sub_matrices()
{
  int i, nnz, *rowbeg, *colidx;
  double *vals;
  //
  // reform dofR, dofC, and refactor AR if extra corners were added in
  // determine_extra_corners
  //
  
  nR = 0; nC = 0;
  for (i=0; i<ndof_sub; i++) {
    if (corner_flag[i] == 0) nR++;
    else nC++;
  }
  dofR = new int[nR]; dofC = new int[nC];
  nR = 0; nC = 0;
  for (i=0; i<ndof_sub; i++) {
    if (corner_flag[i] == 0) {
      dofR[nR] = i;
      nR++;
    }
    else {
      dofC[nC] = i;
      nC++;
    }
  }
  gen_matrix(Block_Stiff, nR, dofR, rowbeg, colidx, vals);
  AR = new CLAPS_sparse_lu();
  nnz = rowbeg[nR];
  AR->factor(nR, nnz, rowbeg, colidx, vals, scale_option);
  delete [] rowbeg; delete [] colidx; delete [] vals;
  //
  // determine internal and boundary dofs and factor AI
  //
  nI = 0; nB = 0;
  for (i=0; i<ndof_sub; i++) {
    if ((sub2[i+1]-sub2[i]) == 1) nI++;
    else nB++;
  }
  dofI = new int[nI]; dofB = new int[nB];
  nI = 0; nB = 0;
  for (i=0; i<ndof_sub; i++) {
    if ((sub2[i+1]-sub2[i]) == 1) {
      dofI[nI] = i;
      nI++;
    }
    else {
      dofB[nB] = i;
      nB++;
    }
  }
  assert ((nR + nC) == (nI + nB));
  gen_matrix(Block_Stiff, nI, dofI, rowbeg, colidx, vals);
  AI = new CLAPS_sparse_lu();
  nnz = rowbeg[nI];
  AI->factor(nI, nnz, rowbeg, colidx, vals, scale_option);
  delete [] rowbeg; delete [] colidx; delete [] vals;
}

void CLIP_solver::calculate_coarse()
{
  int i, j, ibeg, dof, *sub_to_R, INFO, NumEntries, *Indices, gmin, base;
  double *sub_diag, *CT, sum, *Values;

  sub_to_R = new int[ndof_sub]; for (i=0; i<ndof_sub; i++) sub_to_R[i] = -1;
  for (i=0; i<nR; i++) sub_to_R[dofR[i]] = i;

  //
  // specify global numbers for columns of Phi
  //
  int *columns = new int[nC+ndof_set];
  for (i=0; i<nC; i++) columns[i] = SubMap->GID(dofC[i]);
  base = StandardMap->MaxAllGID() + 1;
  for (i=0; i<ndof_set; i++) {
    gmin = base;
    for (j=dset2[i]; j<dset2[i+1]; j++) {
      dof = dset1[j];
      if (SubMap->GID(dof) < gmin) gmin = SubMap->GID(dof);
    }
    columns[nC+i] = gmin;
  }
  Epetra_Map Columns(-1, nC+ndof_set, columns, 0, Comm);
  //
  // specify global numbers of columns of Phi owned by processor
  //
  int ncol_own(0);
  for (i=0; i<(nC+ndof_set); i++) {
    dof = SubMap->LID(columns[i]);
    if (owner_flag[dof] == true) ncol_own++;
  }
  int *columns_own = new int[ncol_own];
  ncol_own = 0;
  for (i=0; i<(nC+ndof_set); i++) {
    dof = SubMap->LID(columns[i]);
    if (owner_flag[dof] == true) {
      columns_own[ncol_own] = columns[i];
      ncol_own++;
    }
  }
  Epetra_Map Columns_own(-1, ncol_own, columns_own, 0, Comm);
  delete [] columns; delete [] columns_own;
  //
  // determine entries of constraint matrix and substructure dof weights
  //
  get_matrix_diag(Block_Stiff, sub_diag);
  Epetra_Vector Sub_diag(View, *SubMap, sub_diag);
  Epetra_Vector Glob_diag(*StandardMap);
  Epetra_Vector Sub_diag_sum(*SubMap);
  Glob_diag.Export(Sub_diag, *Exporter, Add);
  Sub_diag_sum.Import(Glob_diag, *Importer, Insert);
  weight = new double[nB];
  for (i=0; i<nB; i++) weight[i] = sub_diag[dofB[i]]/Sub_diag_sum[dofB[i]];
  CT = new double[nR*ndof_set];
  myzero(CT, nR*ndof_set);
  for (i=0; i<ndof_set; i++) {
    ibeg = i*nR;
    sum = 0;
    for (j=dset2[i]; j<dset2[i+1]; j++) {
      dof = sub_to_R[dset1[j]];
      assert(dof >= 0);
      CT[ibeg+dof] = Sub_diag_sum[dset1[j]];
      sum += CT[ibeg+dof];
    }
    for (j=dset2[i]; j<dset2[i+1]; j++) {
      dof = sub_to_R[dset1[j]];
      CT[ibeg+dof] /= sum;
    }
  }
  delete [] sub_diag;
  /*  
  if (MyPID == 3) {
    std::cout << "weights = " << std::endl;
    for (i=0; i<nB; i++) std::cout << weight[i] << std::endl;
  }
  */
  //
  // solve for substructure coarse interpolation matrices
  //
  ARinvCT = new double[nR*ndof_set];
  double *TEMP = new double[nR*ndof_set];
  if (ndof_set > 0) AR->sol(ndof_set, CT, ARinvCT, TEMP);
  CARinvCT = new double[ndof_set*ndof_set];
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  char TRANSA = 'T'; char TRANSB = 'N'; double ALPHA(1); double BETA(0);
  char UPLO = 'U';
  if (ndof_set > 0)
    EB.GEMM(TRANSA, TRANSB, ndof_set, ndof_set, nR, ALPHA, CT, nR,
	    ARinvCT, nR, BETA, CARinvCT, ndof_set);
  /*
  if (MyPID == 0) {
    std::cout << "CARinvCT = " << std::endl;
    for (i=0; i<ndof_set; i++) {
      for (j=0; j<ndof_set; j++) std::cout << CARinvCT[j*ndof_set+i] << " ";
      std::cout << std::endl;
    }
  }
  */
  if (ndof_set > 0) {
    EL.POTRF(UPLO, ndof_set, CARinvCT, ndof_set, &INFO);
    assert (INFO == 0);
  }
  delete [] TEMP;
  RHS_cg  = new double[ndof_sub]; 
  SOL_cg  = new double[ndof_sub]; 
  TEMP_cg = new double[ndof_sub];
  lambda_e = new double[ndof_set];
  double *phi = new double[ndof_sub];
  int *colPhiT = new int[ndof_sub]; for (i=0; i<ndof_sub; i++) colPhiT[i] = i;
  Epetra_CrsMatrix PhiTs(Copy, Columns, Block_Stiff->RowMap(), ndof_sub, true);
  for (i=0; i<(nC+ndof_set); i++) {
    myzero(phi, ndof_sub);
    if (i < nC) {
      phi[dofC[i]] = 1.0;
      myzero(RHS_cg, nR);
      Block_Stiff->ExtractMyRowView(dofC[i], NumEntries, Values, Indices);
      for (j=0; j<NumEntries; j++) {
	dof = sub_to_R[Indices[j]];
	if (dof >= 0) RHS_cg[dof] = Values[j];
      }
      if (nR > 0) AR->sol(1, RHS_cg, SOL_cg, TEMP_cg);
      if (ndof_set > 0) 
	EB.GEMV(TRANSA, nR, ndof_set, ALPHA, CT, nR, SOL_cg, BETA, lambda_e);
    }
    else {
      myzero(lambda_e, ndof_set);
      lambda_e[i-nC] = 1;
      myzero(SOL_cg, nR);
    }
    if (ndof_set > 0) {
      EL.POTRS(UPLO, ndof_set, 1, CARinvCT, ndof_set, lambda_e, ndof_set, 
	       &INFO);
      assert (INFO == 0);
      EB.GEMV(TRANSB, nR, ndof_set, ALPHA, ARinvCT, nR, lambda_e, -1.0, 
	      SOL_cg);
    }
    for (j=0; j<nR; j++) phi[dofR[j]] = SOL_cg[j];
    PhiTs.InsertMyValues(i, ndof_sub, phi, colPhiT);
  }
  delete [] CT; delete [] sub_to_R; delete [] phi; delete [] colPhiT;
  PhiTs.FillComplete(Block_Stiff->RowMap(), Columns_own);
  PhiTs.OptimizeStorage();
  Epetra_CrsMatrix *PhiT = new Epetra_CrsMatrix(Copy, Columns_own, 0);
  Epetra_Export Exporter_temp(Columns, Columns_own);
  PhiT->Export(PhiTs, Exporter_temp, Insert);
  PhiT->FillComplete(Block_Stiff->RowMap(), Columns_own);
  PhiT->OptimizeStorage();
  Epetra_CrsMatrix *KPhi, *Phi, *Kc_dist;
  construct_transpose(PhiT, Phi);
  Phi->OptimizeStorage();
  EpetraExtCD::MatrixMatrix::Multiply(*Block_Stiff, *Phi, KPhi);
  EpetraExtCD::MatrixMatrix::Multiply(*PhiT, *KPhi, Kc_dist);
  delete PhiT; delete KPhi;
  //
  // reduce Phi to boundary dofs
  //
  int *gdofB = new int[nB];
  for (i=0; i<nB; i++) gdofB[i] = Block_Stiff->GRID(dofB[i]);
  Epetra_Map Map_temp(-1, nB, gdofB, 0, Comm);
  delete [] gdofB;
  PhiB = new Epetra_CrsMatrix(Copy, Map_temp, 0);
  Epetra_Import Importer_temp(Map_temp, Block_Stiff->RowMap());
  PhiB->Import(*Phi, Importer_temp, Insert);
  PhiB->FillComplete(Columns_own, Map_temp);
  PhiB->OptimizeStorage();
  delete Phi;
  //  std::cout << *PhiB << std::endl;
  /*  
  Epetra_Vector xB(Columns_own); xB.PutScalar(1.0);
  Epetra_Vector PhiBxB(MapB);
  PhiB->Multiply(false, xB, PhiBxB);
  std::cout << PhiBxB << std::endl;
  */
  //
  // gather coarse stiffness matrix to processor 0
  //
  Epetra_Map Map_Sub(-1, ncol_own, 0, Comm);
  Epetra_IntVector Int_Sub(Map_Sub);
  for (i=0; i<ncol_own; i++) Int_Sub[i] = Columns_own.GID(i);
  int ngather(0);
  if (MyPID == 0) ngather = Map_Sub.NumGlobalPoints();
  Epetra_Map Map_Gather(-1, ngather, 0, Comm);
  Epetra_IntVector Int_Gather(Map_Gather);
  Epetra_Import Import_Gather(Map_Gather, Map_Sub);
  Int_Gather.Import(Int_Sub, Import_Gather, Insert);
  int *gcdof = new int[ngather];
  for (i=0; i<ngather; i++) gcdof[i] = Int_Gather[i];
  Epetra_Map Map_Final(-1, ngather, gcdof, 0, Comm);
  delete [] gcdof;
  Importer_Kc = new Epetra_Import (Map_Final, Kc_dist->RowMap());
  Exporter_Kc = new Epetra_Export (Map_Final, Kc_dist->RowMap());
  Epetra_CrsMatrix Kc(Copy, Map_Final, 0);
  Kc.Import(*Kc_dist, *Importer_Kc, Insert);
  Kc.FillComplete();
  Kc.OptimizeStorage();
  work_Kc = new Epetra_Vector(Kc.RowMap());
  delete Kc_dist;
  //  std::cout << Kc << std::endl;
  /*    
  Epetra_Vector X(Kc.RowMap());
  X.PutScalar(1.0);
  Epetra_Vector Kc_X(Kc.RowMap());
  Kc.Multiply(false, X, Kc_X);
  std::cout << Kc_X << std::endl;
  */
  double *Xvecs = 0;
  //
  // factor coarse stiffness matrix
  //
  if (MyPID == 0) {
    int *rowbeg, *colidx, nnz;
    double *vals;
    rowbeg = new int[ngather+1]; rowbeg[0] = 0;
    Kc.ExtractMyRowView(0, NumEntries, vals, colidx);
    for (i=0; i<ngather; i++) {
      Kc.ExtractMyRowView(i, NumEntries, Values, Indices);
      rowbeg[i+1] = rowbeg[i] + NumEntries;
    }
    if (print_flag > 0) 
      fout << "coarse problem dimension = " << ngather << std::endl;
    nnz = rowbeg[ngather];
    /*
    std::ofstream ffout;
    ffout.open("coarse_mat.dat");
    for (i=0; i<ngather; i++) 
      for (j=rowbeg[i]; j<rowbeg[i+1]; j++)
	ffout << i+1 << " " << colidx[j]+1 << " " << vals[j] << std::endl;
    ffout.close();
    */
    CLAPS_sparse_lu *AA;
    num_tied_down = 0;
    if (num_rigid_mode > ngather) num_rigid_mode = ngather;
    if (num_rigid_mode > 0) {
      CRD_utils::tie_down_coarse(ngather, rowbeg, colidx, vals,
				 num_rigid_mode, scale_option, num_tied_down, 
				 tied_down, AA, Xvecs);
    }
    if (print_flag > 0) fout << "num_rigid_mode, num_tied_down = " 
			     << num_rigid_mode << " "
			     << num_tied_down << std::endl;
    AKc = new CLAPS_sparse_lu();
    AKc->factor(ngather, nnz, rowbeg, colidx, vals, scale_option);
    delete [] rowbeg;
    SOL_Kc  = new double[ngather];
    TEMP_Kc = new double[ngather];
  }
  //
  // set up maps, importers, and exporters
  //
  nB_own = 0;
  gdofB = new int[nB];
  for (i=0; i<nB; i++) {
    gdofB[i] = SubMap->GID(dofB[i]);
    if (owner_flag[dofB[i]] == true) nB_own++;
  }
  int *gdofB_own = new int[nB_own];
  nB_own = 0;
  for (i=0; i<nB; i++) {
    if (owner_flag[dofB[i]] == true) {
      gdofB_own[nB_own] = SubMap->GID(dofB[i]);
      nB_own++;
    }
  }
  MapB_Sub = new Epetra_Map(-1, nB, gdofB, 0, Comm);
  MapB_St = new Epetra_Map(-1, nB_own, gdofB_own, 0, Comm);
  delete [] gdofB; delete [] gdofB_own;
  ImporterB = new Epetra_Import(*MapB_Sub, *MapB_St);
  ExporterB = new Epetra_Export(*MapB_Sub, *MapB_St);
  //
  // calculate rigid body modes for mode acceleration technique
  //
  Comm.Broadcast(&num_tied_down, 1, 0);  
  if (num_tied_down > 0) {
    Epetra_MultiVector Evec(View, Map_Final, Xvecs, ngather, num_tied_down);
    Epetra_MultiVector Evec_B(PhiB->DomainMap(), num_tied_down);
    Evec_B.Export(Evec, *Exporter_Kc, Insert);
    Epetra_MultiVector Phir_Sub(PhiB->RowMap(), num_tied_down);
    Epetra_MultiVector Phir_Sub2(*MapB_Sub, num_tied_down);
    PhiB->Multiply(false, Evec_B, Phir_Sub);
    int MyLDA, MyLDA2, jbeg, jbeg2;
    double *phirsub, *phirsub2;
    Phir_Sub.ExtractView(&phirsub, &MyLDA);
    Phir_Sub2.ExtractView(&phirsub2, &MyLDA2);
    for (j=0; j<num_tied_down; j++) {
      jbeg = j*MyLDA;
      jbeg2 = j*MyLDA2;
      for (i=0; i<nB; i++) phirsub2[jbeg2+i] = phirsub[jbeg+i]*weight[i];
    }
    Phir_St = new Epetra_MultiVector(*MapB_St, num_tied_down);
    Phir_St->Export(Phir_Sub2, *ExporterB, Add);

    double *infnorm = new double[num_tied_down];
    Phir_St->NormInf(infnorm);
    double *phirst;
    Phir_St->ExtractView(&phirst, &MyLDA);
    for (i=0; i<num_tied_down; i++)
      for (j=0; j<Phir_St->MyLength(); j++) phirst[j+i*MyLDA] /= infnorm[i];
    delete [] infnorm;

    Epetra_LocalMap Local_Map(num_tied_down, 0, Comm);
    Epetra_MultiVector AA(Local_Map, num_tied_down);
    Rhs_null = new Epetra_Vector(Local_Map);
    AA.Multiply('T', 'N', 1.0, *Phir_St, *Phir_St, 0.0);
    double *aa;
    AA.ExtractView(&aa, &MyLDA);
    PhirTPhir = new double[num_tied_down*num_tied_down];
    for (i=0; i<num_tied_down; i++)
      for (j=0; j<num_tied_down; j++) 
	PhirTPhir[i+j*num_tied_down] = aa[i+j*MyLDA];
    EL.POTRF(UPLO, num_tied_down, PhirTPhir, num_tied_down, &INFO);
    assert (INFO == 0);
  }
  delete [] Xvecs;
}

void CLIP_solver::get_matrix_diag(Epetra_CrsMatrix* A, double* & diag)
{
  int i, j, NumEntries, *Indices, n, dof;
  double *Values;
  n = A->NumMyRows();
  diag = new double[n]; myzero(diag, n);
  for (i=0; i<n; i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) if (Indices[j] == i) diag[i] = Values[j];
  }
}

void CLIP_solver::gen_matrix(Epetra_CrsMatrix* A, int na, int adof[], 
			     int* & rowbeg, int* & colidx, double* &vals)
{
  int i, j, NumEntries, *Indices, n, *imap, nnz(0), dof;
  double *Values;
  n = A->NumMyRows();
  imap = new int[n]; for (i=0; i<n; i++) imap[i] = -1;
  rowbeg = new int[na+1]; rowbeg[0] = 0;
  for (i=0; i<na; i++) imap[adof[i]] = i;
  for (i=0; i<na; i++) {
    A->ExtractMyRowView(adof[i], NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) if (imap[Indices[j]] >= 0) nnz++;
    rowbeg[i+1] = nnz;
  }
  colidx = new int[nnz];
  vals = new double[nnz];
  nnz = 0;
  for (i=0; i<na; i++) {
    A->ExtractMyRowView(adof[i], NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      dof = Indices[j];
      if (imap[dof] >= 0) {
	colidx[nnz] = imap[dof];
	vals[nnz] = Values[j];
	nnz++;
      }
    }
  }
  delete [] imap; 
}

int CLIP_solver::determine_shell(int ii)
{
  int i, dof, flag(0), *loc_dofs;
  NodalDofs->ExtractView(&loc_dofs);
  for (i=dset2[ii]; i<dset2[ii+1]; i++) {
    dof = dset1[i];
    if ((loc_dofs[dof] >= 4) && (loc_dofs[dof] <= 6)) flag = 1;
  }
  return flag;
}

void CLIP_solver::get_two_points(double xyz[], int LDA, int adof[],
				 int ndof_cand, int &i1, int &i2)
{
  int i, j, dof1, dof2;
  double dist2, max_dist2(-1), x1, y1, z1, x2, y2, z2;
  for (i=0; i<ndof_cand; i++) {
    dof1 = adof[i];
    x1 = xyz[dof1]; y1 = xyz[dof1+LDA]; z1 = xyz[dof1+2*LDA];
    for (j=i+1; j<ndof_cand; j++) {
      dof2 = adof[j];
      x2 = xyz[dof2]; y2 = xyz[dof2+LDA]; z2 = xyz[dof2+2*LDA];
      dist2 = calc_dist2(x1, y1, z1, x2, y2, z2);
      if (dist2 > max_dist2) {
	max_dist2 = dist2;
	i1 = dof1;
	i2 = dof2;
      }
    }
  }
  assert (max_dist2 >= 0);
}

void CLIP_solver::get_one_point(double xyz[], int LDA, int adof[],
				int ndof_cand, int i1, int &i2)
{
  int i, dof;
  double dist2, max_dist2(-1), x1, y1, z1, x2, y2, z2;
  x1 = xyz[i1]; y1 = xyz[i1+LDA]; z1 = xyz[i1+2*LDA];
  for (i=0; i<ndof_cand; i++) {
    dof = adof[i];
    if (dof != i1) {
      x2 = xyz[dof]; y2 = xyz[dof+LDA]; z2 = xyz[dof+2*LDA];
      dist2 = calc_dist2(x1, y1, z1, x2, y2, z2);
      if (dist2 > max_dist2) {
	max_dist2 = dist2;
	i2 = dof;
      }
    }
  }
  assert (max_dist2 >= 0);
}

void CLIP_solver::get_cand_dofs(int ii, int adof[], bool adof_flag[], 
				int &ndof_cand)
{
  int i, j, k, dof, subnum, loc_dof, loc_dof1, *loc_dofs, nsuba, flag;
  NodalDofs->ExtractView(&loc_dofs);
  dof = dset1[dset2[ii]];
  loc_dof1 = loc_dofs[dof];
  nsuba = sub2[dof+1] - sub2[dof];
  int *suba = new int[nsuba];
  for (i=sub2[dof]; i<sub2[dof+1]; i++) suba[i-sub2[dof]] = sub1[i];
  int *flaga = new int[nsuba];
  ndof_cand = 0;
  for (i=0; i<ndof_set; i++) {
    dof = dset1[dset2[i]];
    for (j=0; j<nsuba; j++) flaga[j] = 0;
    if ((sub2[dof+1]-sub2[dof]) > 1) {
      for (j=sub2[dof]; j<sub2[dof+1]; j++) {
	subnum = sub1[j];
	for (k=0; k<nsuba; k++) {
	  if (subnum == suba[k]) flaga[k] = 1;
	}
      }
      flag = 1;
      for (j=0; j<nsuba; j++) if (flaga[j] == 0) flag = 0;
      if (flag == 1) {
	for (j=dset2[i]; j<dset2[i+1]; j++) {
	  dof = dset1[j];
	  loc_dof = loc_dofs[dof];
	  if ((adof_flag[dof] == false) && (loc_dof == loc_dof1)) {
	    adof_flag[dof] = true;
	    adof[ndof_cand] = dof;
	    ndof_cand++;
	  }
	}
      }
    }
  }
  for (i=0; i<ndof_cand; i++) adof_flag[adof[i]] = false;
  delete [] suba; delete [] flaga;
}

void CLIP_solver::get_first_point(int adof[], int ndof_cand, int &i1)
{
  int i, dof, *gdofs, max_gdof(-7);
  gdofs = SubMap->MyGlobalElements();
  for (i=0; i<ndof_cand; i++) {
    dof = adof[i];
    if (gdofs[dof] > max_gdof) {
      i1 = dof;
      max_gdof = gdofs[dof];
    }
  }
  return;
}

void CLIP_solver::get_third_point(double xyz[], int LDA, int adof[],
				  int ndof_cand, int i1, int i2, int &i3)
{
  int i, dof;
  double dx1, dy1, dz1, dx2, dy2, dz2, area2, max_area2(-7);
  dx1 = xyz[i2      ] - xyz[i1      ];
  dy1 = xyz[i2+LDA  ] - xyz[i1+LDA  ];
  dz1 = xyz[i2+2*LDA] - xyz[i1+2*LDA];
  for (i=0; i<ndof_cand; i++) {
    dof = adof[i];
    dx2 = xyz[dof      ] - xyz[i1      ];
    dy2 = xyz[dof+LDA  ] - xyz[i1+LDA  ];
    dz2 = xyz[dof+2*LDA] - xyz[i1+2*LDA];
    area2 = calc_area2(dx1, dy1, dz1, dx2, dy2, dz2);
    if (area2 > max_area2) {
      i3 = dof;
      max_area2 = area2;
    }
  }
}

double CLIP_solver::calc_dist2(double x1, double y1, double z1,
			       double x2, double y2, double z2)
{
  double dx, dy, dz, dist2;
  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  dist2 = dx*dx + dy*dy + dz*dz;
  return dist2;
}

double CLIP_solver::calc_area2(double dx1, double dy1, double dz1,
			       double dx2, double dy2, double dz2)
{
  double dx3, dy3, dz3, area2;
  dx3 = dy1*dz2 - dz1*dy2;
  dy3 = dz1*dx2 - dx1*dz2;
  dz3 = dx1*dy2 - dy1*dx2;
  area2 = dx3*dx3 + dy3*dy3 + dz3*dz3;
  return area2;
}

void CLIP_solver::determine_components(int A1[], int A2[], int N, 
				       int* & compt1, int* & compt2, int & ncompt)
{
  int i, ic;
  CRD_utils::Graph_class Graph(N, A1, A2);
  int *component = new int[N];
  Graph.Components(component, ncompt);
  compt1 = new int[N];
  compt2 = new int[ncompt+1]; compt2[0] = 0;
  int *count_comp = new int[ncompt]; myzero(count_comp, ncompt);
  for (i=0; i<N; i++) count_comp[component[i]]++;
  for (i=0; i<ncompt; i++) compt2[i+1] = compt2[i] + count_comp[i];
  myzero(count_comp, ncompt);
  for (i=0; i<N; i++) {
    ic = component[i];
    compt1[compt2[ic] + count_comp[ic]] = i;
    count_comp[ic]++;
  }
  delete [] component; delete [] count_comp;
}



void CLIP_solver::zero_pointers()
{
  Importer = 0; ImporterB = 0; ImporterSt2B = 0; Importer_Kc = 0; Exporter = 0;
  ExporterB = 0; Exporter_lam = 0; Exporter_Kc = 0; ImporterStand = 0;
  ExporterStand = 0; Lambda = 0; Lambda_local = 0; ConError = 0; ApB_Sub = 0;
  ApB_St = 0; pB_Sub = 0; pB_St = 0; zB_Sub = 0; zB_St = 0; uB_St = 0; 
  u_St = 0; vStand = 0; r_St = 0; r_Sub = 0; rB_St = 0; rB_Sub = 0;
  work_Sub = 0; work_Kc = 0; rBSub = 0; workB_St = 0; work_St = 0; 
  work_St2 = 0; vec_PhiB = 0; prod_PhiB = 0; vec1_ASub = 0; vec2_ASub = 0;
  wStand = 0; StandardMap = 0; SubMap = 0; RowMapMyCon = 0; MapB_Sub = 0;
  MapB_St = 0; CtT = 0; Tran = 0; Block_Stiff = 0; PhiB = 0; Coords_red = 0;
  NodalDofs_red = 0; AR = 0; AI = 0; AKc = 0; comp1 = 0; comp2 = 0; sub1 = 0;
  sub2 = 0; dset1 = 0; dset2 = 0; corner_flag = 0; mycdof = 0; dofI = 0;
  dofB = 0; dofC = 0; dofR = 0; lambda = 0; lambda_local = 0; weight = 0;
  ARinvCT = 0; CARinvCT = 0; lambda_e = 0; RHS_cg = 0; SOL_cg = 0; TEMP_cg = 0;
  SOL_Kc = 0; TEMP_Kc = 0; rcurra = 0; rhoa = 0; betaa = 0; pApa = 0; Dtri = 0;
  Etri = 0; econa = 0; P_ortho = 0; AP_ortho = 0; orth1 = 0; orth2 = 0; 
  owner_flag = 0; VV = 0; RR = 0; HH = 0; zz = 0; cc = 0; ss = 0; norms = 0; 
  gmres_vec = 0; gmres_sum = 0; tied_down = 0; Phir_St = 0; PhirTPhir = 0;
  Rhs_null = 0;
}

void CLIP_solver::construct_subdomains()
{
  if (print_flag > 9) fout << "in construct_subdomains" << std::endl;
}

int CLIP_solver::initialize_subdomains()
{
  if (print_flag > 9) fout << "in initialize_subdomains" << std::endl;
  return(0);
}

void CLIP_solver::calculate_coarse_stiff()
{
}

void CLIP_solver::gather_coarse_stiff()
{
}

void CLIP_solver::factor_coarse_stiff()
{
}

void CLIP_solver::solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
			int & num_iter, int & solver_status,
			int & max_added_cor)
{
  int pcg_status, gmres_status;
  double starttime, endtime;
  max_added_cor = max_added_corner;
  Comm.Barrier();
  if (MyPID == 0) starttime = MPI_Wtime();
  if (print_flag >= 0) fout.open("CLIP_solver.data", std::ios::app);
  if (krylov_method != 1) {
    pcg_solve(uStand, fStand, num_iter, pcg_status);
    solver_status = pcg_status;
  }
  if (krylov_method == 1) {
    gmres_solve(uStand, fStand, num_iter, gmres_status);
    solver_status = gmres_status;
  }
  Comm.Barrier();
  if (print_flag > 0) {
    endtime = MPI_Wtime();
    fout << "elapsed time for clip solver solve = " << endtime-starttime
	 << " seconds" << std::endl;
  }
  if (print_flag >= 0) fout.close();
}

void CLIP_solver::pcg_init()
{
  pB_Sub  = new Epetra_Vector(*MapB_Sub);
  pB_St   = new Epetra_Vector(*MapB_St);
  zB_Sub  = new Epetra_Vector(*MapB_Sub);
  zB_St   = new Epetra_Vector(*MapB_St);
  ApB_Sub = new Epetra_Vector(*MapB_Sub);
  ApB_St  = new Epetra_Vector(*MapB_St);
  uB_St   = new Epetra_Vector(*MapB_St);
  ImporterSt2B = new Epetra_Import(*MapB_St, *StandardMap);

  r_St     = new Epetra_Vector(*StandardMap);
  rB_St    = new Epetra_Vector(*MapB_St);
  r_Sub    = new Epetra_Vector(*SubMap);
  rB_Sub   = new Epetra_Vector(*MapB_Sub);

  work_Sub  = new Epetra_Vector(*SubMap);
  work_St   = new Epetra_Vector(*StandardMap);
  work_St2  = new Epetra_Vector(*StandardMap);
  workB_St  = new Epetra_Vector(*MapB_St);
  vec_PhiB  = new Epetra_Vector(PhiB->RowMap());
  prod_PhiB = new Epetra_Vector(PhiB->DomainMap());

  rcurra = new double[maxiter+2];

  int nrB_St = rB_St->GlobalLength();
  if (nrB_St < max_orthog) max_orthog = nrB_St; 

  if (krylov_method != 1) {
    rhoa = new double[maxiter+2];
    betaa = new double[maxiter+2];
    pApa = new double[maxiter+2];
    Dtri = new double[maxiter+2];
    Etri = new double[maxiter+2];
    econa = new double[maxiter+2];
    P_ortho  = new double[max_orthog*nB_own];
    AP_ortho = new double[max_orthog*nB_own];
    orth1 = new double[max_orthog];
    orth2 = new double[max_orthog];
  }
  if (krylov_method == 1) {
    int m = nB_own*(maxiter+1);
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
  cg_iter = -1;
  n_orthog_used = 0;
}

void CLIP_solver::pcg_solve(Epetra_Vector* uStand, const Epetra_Vector* fStand,
			    int & num_iter, int & pcg_status)
{
  int i, iflag(1), max_iter_new;
  double dprod, rorig, beta, roldzold(0), pAp, alpha, rcurr, ractual;
  double norm_rconstraint, norm_conerror, tpe, rnorm;
  //
  // determine residual for constrained problem
  //
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *r_St);
  if (ncon_global == 0) r_St->Update(1.0, *fStand, 0.0);
  uB_St->PutScalar(0);
  //
  // calculate initial residual
  //
  r_St->Norm2(&rorig);
  rcurra[0] = rorig;
  if (print_flag > 0) {
    fout << "original residual                          = " << rorig << std::endl;
  }
  if (rorig == 0) {
    uStand->PutScalar(0);
    num_iter = 0;
    pcg_status = 0;
    return;
  }
  /*
  if (num_tied_down > 0) {
    int *nodaldofs;
    double *rst, sum, sum_all, rn;
    r_St->ExtractView(&rst);
    Epetra_IntVector NDOF(*StandardMap);
    NDOF.Export(*NodalDofs, *Exporter, Insert);
    NDOF.ExtractView(&nodaldofs);
    sum = 0;
    assert (r_St->MyLength() == NDOF.MyLength()); 
    for (i=0; i<r_St->MyLength(); i++)
      if (nodaldofs[i] == 1) sum += rst[i];
    Comm.SumAll(&sum, &sum_all, 1);
    r_St->Norm2(&rn);
    sum_all /= rn;
    if (MyPID == 0) std::cout << "1-direction normalized force = " << sum_all
			 << std::endl;
  }
  */
  //
  // calculate residual at boundary and total potential energy following 
  // initial static condensation
  //
  tpe = stat_cond();
  rB_St->Norm2(&rnorm);
  if (print_flag > 0) {
    fout << "----------after initial static condensation---------" << std::endl;
    fout << "residual                                   = " << rnorm << std::endl;
    fout << "total potential energy reduction           = " << tpe << std::endl;
    if (n_orthog_used == 0)
      fout << "----------------------------------------------------" << std::endl;
  }
  //
  // remove part of residual parallel to null space for mode acceleration
  // technique
  //
  if (num_tied_down > 0) {
    remove_orthog_null(rB_St);
    if (print_flag > 0) fout << "singular system being solved" << std::endl;
  }
  //
  // reduce total potential energy further by using stored vectors
  //
  pcg_status = 1;
  max_iter_new = maxiter;
  if (n_orthog_used > 0) {
    tpe = initial_update();
    rB_St->Norm2(&rnorm);
    if (print_flag > 0) {
      fout << "----------after using stored search directions------" << std::endl;
      fout << "number of search directions used           = " 
	   << n_orthog_used << std::endl;
      fout << "residual                                   = " << rnorm
	   << std::endl;
      fout << "additional reduction in tpe                = " << tpe << std::endl;
      if (tpe > 0) {
	fout << "Warning: using stored search directions increased ";
	fout << "         total potential energy" << std::endl;
      }
      fout << "----------------------------------------------------" << std::endl;
    }
    if (rnorm/rorig <= solver_tol) {
      pcg_status = max_iter_new = num_iter = 0;
    }
  }
  //
  // pcg iterations
  //
  if (n_orthog_used == max_orthog) cg_iter = 0;
  for (int iter=0; iter<max_iter_new; iter++) {
    if (MyPID == -1) std::cout << "iteration " << iter+1 << " of maxiter = "
    			 << maxiter << std::endl;
    //
    // apply BDDC preconditioner
    //
    apply_preconditioner();
    //
    // augment standard BDDC preconditioner with stored vectors if in
    // standard pcg mode
    //
    if ((cg_iter >= 0) && (n_orthog_used > 0)) {
      determine_AxB(zB_St, ApB_St);
      ApB_St->Update(1.0, *rB_St, -1.0);
      final_update(zB_St, ApB_St);
    }
    rB_St->Dot(*zB_St, &dprod);
    if (krylov_method == 0) assert (dprod >= 0);
    if (cg_iter >= 0) {
      rhoa[cg_iter] = sqrt(fabs(dprod));
      if (cg_iter == 0) {
	beta = 0;
	pB_St->Update(1.0, *zB_St, 0.0);
      }
      else {
	beta = dprod/roldzold;
	pB_St->Update(1.0, *zB_St, beta);
      }
      betaa[cg_iter] = beta;
      roldzold = dprod;
    }
    else {
      pB_St->Update(1.0, *zB_St, 0.0);
    }
    //
    // project off stored vectors from preconditioned residual
    //
    if ((n_orthog_used > 0) && (cg_iter == -1)) remove_search(pB_St);
    //
    // determine ABB * pB
    //
    determine_AxB(pB_St, ApB_St);
    pB_St->Dot(*ApB_St, &pAp);
    if (cg_iter >= 0) {
      pApa[cg_iter] = pAp;
      cg_iter++;
    }
    alpha = dprod/pAp;
    //
    // store pB_St and ApB_St
    //
    store_search(pAp);
    //
    // update uB_St and rB_St
    //
    uB_St->Update( alpha, *pB_St,  1.0);
    rB_St->Update(-alpha, *ApB_St, 1.0);
    rB_St->Dot(*rB_St, &rcurr);
    rcurr = sqrt(rcurr);
    rcurra[iter+1] = rcurr;
    num_iter = iter+1;
    if (MyPID == -1) std::cout << "n_orthog_used, cg_iter = " << n_orthog_used
			 << " " << cg_iter << std::endl;
    if (MyPID == -1) std::cout << "alpha, dprod, rtol = " << alpha << " " 
    			 << dprod << " " << rcurr/rorig << std::endl;
    if ((iflag > 0) && (rcurr/rorig <= solver_tol)) {
      pcg_status = 0;
      break;
    }
  }
  //
  // calculate uStand
  //
  if (ncon_global == 0) r_St->Update(1.0, *fStand, 0.0);
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *r_St);
  calculate_u_St();
  if (ncon_global == 0) uStand->Update(1.0, *u_St, 0.0);
  if (ncon_global > 0) Tran->Multiply(false, *u_St, *uStand);
  //
  // calculate actual residual for reduced system
  //
  calculate_Au(u_St, work_St2);
  work_St2->Update(1.0, *r_St, -1.0);
  work_St2->Dot(*work_St2, &ractual);
  ractual = sqrt(ractual);
  if ((ractual/rorig > 10*solver_tol) && (num_tied_down == 0)) pcg_status = 1;
  //
  // calculate Lagrange multipliers and check residual for full system
  //
  if (ncon_global > 0) {
    calculate_AuStand(uStand, vStand);
    vStand->Update(1.0, *fStand, -1.0);
    calculate_multipliers(uStand, norm_rconstraint, norm_conerror);
  }
  if (MyPID == 0) {
    //    std::cout << "rorig                 = " << rorig << std::endl;
    if (print_flag > 0) {
      if (num_iter > 0) fout << "rcurr(recursive)      = " << rcurr << std::endl;
      fout << "rcurr(actual)         = " << ractual << std::endl;
      if (ncon_global > 0) {
	fout << "rcurr(constraint)     = " << norm_rconstraint << std::endl;
	fout << "constraint error norm = " << norm_conerror << std::endl;
      }
      fout << "number of iterations  = " << num_iter << std::endl;
      fout << "solver tolerance      = " << solver_tol << std::endl;
    }
    if (cg_iter > 0) {
      if ((print_flag == 2) || (print_flag == 12)) {
	fout << "condition # estimate      relative residual" 
	     << "   iteration" << std::endl;
	calculate_condition(cg_iter);
	fout << std::setiosflags(std::ios::scientific | std::ios::uppercase);
	for (i=0; i<num_iter; i++) {
	  double ee = 0;
	  if (i >= (num_iter-cg_iter)) ee = econa[i-(num_iter-cg_iter)]; 
	  fout << " " 
	       << std::setw(17) << std::setprecision(10) << ee 
	       << "       " 
	       << std::setw(17) << std::setprecision(10) << rcurra[i+1]/rorig
	       << "        " 
	       << i+1 << std::endl;
	}
      }
      fout << std::resetiosflags(std::ios::scientific);
      fout << std::resetiosflags(std::ios::uppercase);
      fout << std::setprecision(6);
    }
  }
}

void CLIP_solver::gmres_solve(Epetra_Vector* uStand, 
          const Epetra_Vector* fStand, int & num_iter, int & gmres_status)
{
  int i, iflag(1), gmres_iter;
  double dprod, rorig, normb, rcurr, ractual, rnorm;
  double norm_rconstraint, norm_conerror, *vals;
  //
  // determine residual for constrained problem
  //
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *r_St);
  if (ncon_global == 0) r_St->Update(1.0, *fStand, 0.0);
  uB_St->PutScalar(0);
  //
  // calculate initial residual
  //
  r_St->Norm2(&rorig);
  rcurra[0] = rorig;
  if (print_flag > 0) {
    fout << "original residual                          = " << rorig << std::endl;
  }
  if (rorig == 0) {
    uStand->PutScalar(0);
    num_iter = 0;
    gmres_status = 0;
    return;
  }
  //
  // calculate residual at boundary following initial static condensation
  //
  stat_cond();
  rB_St->Norm2(&rnorm);
  if (print_flag > 0) {
    fout << "----------after initial static condensation---------" << std::endl;
    fout << "residual                                   = " << rnorm << std::endl;
    if (n_orthog_used == 0)
      fout << "----------------------------------------------------" << std::endl;
  }
  //
  // remove part of residual parallel to null space for mode acceleration
  // technique
  //
  if (num_tied_down > 0) {
    remove_orthog_null(rB_St);
    if (print_flag > 0) fout << "singular system being solved" << std::endl;
  }
  //
  // gmres iterations
  //
  normb = rnorm;
  gmres_status = 1;
  rB_St->ExtractView(&vals);
  for (i=0; i<nB_own; i++) {
    vals[i] /= normb;
    VV[i] = vals[i];
  }
  for (gmres_iter=0; gmres_iter<maxiter; gmres_iter++) {
    if (MyPID == -1) std::cout << "iteration " << gmres_iter+1 
			  << " of maxiter = " << maxiter << std::endl;
    apply_preconditioner();
    //
    // gmres stuff
    //
    determine_AxB(zB_St, rB_St);
    //
    // two steps of classical Gram-Schmidt
    //
    two_steps_CGS(gmres_iter, rB_St);
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
    memcpy(vals, &VV[nB_own*num_iter], nB_own*sizeof(double));
  }
  //
  // construct the solution
  //
  construct_solution(num_iter, normb);
  //
  // calculate uStand
  //
  if (ncon_global == 0) r_St->Update(1.0, *fStand, 0.0);
  if (ncon_global >  0) Tran->Multiply(true, *fStand, *r_St);
  calculate_u_St();
  if (ncon_global == 0) uStand->Update(1.0, *u_St, 0.0);
  if (ncon_global > 0) Tran->Multiply(false, *u_St, *uStand);
  //
  // calculate actual residual for reduced system
  //
  calculate_Au(u_St, work_St2);
  work_St2->Update(1.0, *r_St, -1.0);
  work_St2->Dot(*work_St2, &ractual);
  ractual = sqrt(ractual);
  if (ractual/rorig > 10*solver_tol) gmres_status = 1;
  //
  // calculate Lagrange multipliers and check residual for full system
  //
  if (ncon_global > 0) {
    calculate_AuStand(uStand, vStand);
    vStand->Update(1.0, *fStand, -1.0);
    calculate_multipliers(uStand, norm_rconstraint, norm_conerror);
  }
  if (MyPID == 0) {
    //    std::cout << "rorig                 = " << rorig << std::endl;
    if (print_flag > 0) {
      if (num_iter > 0) fout << "rcurr(recursive)      = " << rcurr << std::endl;
      fout << "rcurr(actual)         = " << ractual << std::endl;
      if (ncon_global > 0) {
	fout << "rcurr(constraint)     = " << norm_rconstraint << std::endl;
	fout << "constraint error norm = " << norm_conerror << std::endl;
      }
      fout << "number of iterations  = " << num_iter << std::endl;
      fout << "solver tolerance      = " << solver_tol << std::endl;
    }
  }
}

void CLIP_solver::remove_orthog_null(Epetra_Vector *vec)
{
  int i, MyLDA, INFO;
  double *rhsnull, rnorm, rhsnorm, *phirst, *zbst;
  Epetra_LAPACK EL;
  char UPLO = 'U';
  /*
  Phir_St->ExtractView(&phirst, &MyLDA);
  zB_St->ExtractView(&zbst);
  for (i=0; i<nB_own; i++) zbst[i] = phirst[i];
  determine_AxB(zB_St, ApB_St);
  zB_St->Norm2(&rnorm);
  ApB_St->Scale(1/rnorm/Block_Stiff->NormInf());
  std::cout << *ApB_St << std::endl;
  */
  Rhs_null->Multiply('T', 'N', 1.0, *Phir_St, *vec, 0.0);
  Rhs_null->Norm2(&rhsnorm);
  vec->Norm2(&rnorm);
  if (print_flag > 0) fout << "check_orthog_null (before) = " << rhsnorm/rnorm
			   << std::endl;
  Rhs_null->ExtractView(&rhsnull);
  EL.POTRS(UPLO, num_tied_down, 1, PhirTPhir, num_tied_down, rhsnull, 
  	   num_tied_down, &INFO);
  //  vec->Multiply('N', 'N', -1.0, *Phir_St, *Rhs_null, 1.0);
  Rhs_null->Multiply('T', 'N', 1.0, *Phir_St, *vec, 0.0);
  //  std::cout << *Rhs_null << std::endl;
  Rhs_null->Norm2(&rhsnorm);
  vec->Norm2(&rnorm);
  if (print_flag > 0) fout << "check_orthog_null (after ) = " << rhsnorm/rnorm 
			   << std::endl;
}

double CLIP_solver::stat_cond()
{
  int i, j, NumEntries, *Indices, dof;
  double *Values, sum, *rsub, *rbsub, *rbSt, dnorm, tpe(0), tpe_sum;
  r_Sub->ExtractView(&rsub);
  rB_Sub->ExtractView(&rbsub);
  r_Sub->Import(*r_St, *Importer, Insert);
  for (i=0; i<nI; i++) RHS_cg[i] = rsub[dofI[i]];
  if (nI > 0) AI->sol(1, RHS_cg, SOL_cg, TEMP_cg);
  for (i=0; i<nI; i++) tpe -= SOL_cg[i]*RHS_cg[i];
  tpe /= 2;
  Comm.SumAll(&tpe, &tpe_sum, 1);
  myzero(RHS_cg, ndof_sub);
  for (i=0; i<nI; i++) RHS_cg[dofI[i]] = SOL_cg[i];
  for (i=0; i<nB; i++) {
    dof = dofB[i];
    Block_Stiff->ExtractMyRowView(dof, NumEntries, Values, Indices);
    sum = 0;
    for (j=0; j<NumEntries; j++) sum += Values[j]*RHS_cg[Indices[j]];
    rbsub[i] = sum;
  }
  workB_St->PutScalar(0.0);
  workB_St->Export(*rB_Sub, *ExporterB, Add);
  rB_St->Import(*r_St, *ImporterSt2B, Insert);
  rB_St->Update(-1.0, *workB_St, 1.0);
  /* 
  rB_St->Norm2(&dnorm);
  if (MyPID == 0) std::cout << "scale_option = " << scale_option << std::endl;
  if (MyPID == 0) std::cout << "2 norm of rB_St = " << dnorm << std::endl;
  work_Sub->PutScalar(0.0);
  double *worksub;
  work_Sub->ExtractView(&worksub);
  for (i=0; i<ndof_sub; i++) {
    Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
    sum = 0;
    for (j=0; j<NumEntries; j++) sum += Values[j]*RHS_cg[Indices[j]];
    worksub[i] = sum;
  }
  if (MyPID == 0) {
    char fname[101];
    sprintf(fname,"%s%d","test_crd_", MyPID);
    sprintf(fname, "%s.dat", fname);
    std::ofstream ffout;
    ffout.open(fname);
    ffout << ndof_sub << std::endl;
    ffout << nI << std::endl;
    ffout << Block_Stiff->NumMyNonzeros() << std::endl;
    for (i=0; i<ndof_sub; i++) ffout << rsub[i] << std::endl;
    for (i=0; i<nI; i++) ffout << dofI[i] << std::endl;
    for (i=0; i<ndof_sub; i++) {
      Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
      ffout << NumEntries << std::endl;
      for (j=0; j<NumEntries; j++) ffout << Indices[j] << " " 
					 << Values[j] << std::endl;
    }
    ffout.close();
  }
  double *rst, delta;
  r_St->ExtractView(&rst);
  sum = 0;
  for (i=0; i<nI; i++) {
    delta = rst[dofI[i]] - worksub[dofI[i]];
    sum += delta*delta;
  }
  sum = sqrt(sum);
  std::cout << "MyPID, internal norm = " << MyPID << " " << sum << std::endl;

  work_St->PutScalar(0);
  work_St->Export(*work_Sub, *Exporter, Add);
  r_St->Update(-1.0, *work_St, 1.0);
  r_St->Norm2(&dnorm);
  if (MyPID == 0) std::cout << "2 norm of r_St  = " << dnorm << std::endl;
  assert (dnorm == 0);
  */
  return tpe_sum;
}

double CLIP_solver::initial_update()
{
  int i;
  double tpe, *rbst, *pbst, *apbst, dprod1, dprod2;
  Epetra_BLAS EB;
  char TRANS = 'T'; double ALPHA(1); double BETA(0);
  rB_St->ExtractView(&rbst);
  pB_St->ExtractView(&pbst);
  ApB_St->ExtractView(&apbst);
  if ((nB_own > 0) && (n_orthog_used > 0)) 
    EB.GEMV(TRANS, nB_own, n_orthog_used, ALPHA, P_ortho, nB_own, rbst, BETA, 
	    orth1);
  else
    myzero(orth1, n_orthog_used);
  Comm.SumAll(orth1, orth2, n_orthog_used);
  TRANS = 'N';
  if ((nB_own > 0) && (n_orthog_used > 0)) {
    EB.GEMV(TRANS, nB_own, n_orthog_used, ALPHA,  P_ortho, nB_own, orth2, 
	    BETA, pbst);
    EB.GEMV(TRANS, nB_own, n_orthog_used, ALPHA, AP_ortho, nB_own, orth2, 
	    BETA, apbst);
  }
  ApB_St->Dot(*pB_St, &dprod1);
  pB_St->Dot(*rB_St, &dprod2);
  tpe = dprod1/2 - dprod2;
  rB_St->Update(-1.0, *ApB_St, 1.0);
  uB_St->Update( 1.0,  *pB_St, 0.0);
  return tpe;
}

void CLIP_solver::apply_preconditioner()
{
  int i, INFO;
  double *rbsub, *zbsub, *vecphib, *workkc, ALPHA(1.0), BETA(0), nKc;
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  char TRANS = 'T'; char UPLO = 'U';
  zB_Sub->ExtractView(&zbsub);
  rB_Sub->ExtractView(&rbsub);
  vec_PhiB->ExtractView(&vecphib);
  work_Kc->ExtractView(&workkc);
  rB_Sub->Import(*rB_St, *ImporterB, Insert);
  for (i=0; i<nB; i++) rbsub[i] *= weight[i];
  //
  // substructure correction
  //
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nB; i++) TEMP_cg[dofB[i]] = rbsub[i];
  for (i=0; i<nR; i++) RHS_cg[i] = TEMP_cg[dofR[i]];
  if (nR > 0) AR->sol(1, RHS_cg, SOL_cg, TEMP_cg);
  if (ndof_set > 0) {
    EB.GEMV(TRANS, nR, ndof_set, ALPHA, ARinvCT, nR, RHS_cg, BETA, lambda_e);
    EL.POTRS(UPLO, ndof_set, 1, CARinvCT, ndof_set, lambda_e, ndof_set, &INFO);
    assert (INFO == 0);
    ALPHA = -1; BETA = 1; TRANS = 'N';
    EB.GEMV(TRANS, nR, ndof_set, ALPHA, ARinvCT, nR, lambda_e, BETA, SOL_cg);
  }
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nR; i++) TEMP_cg[dofR[i]]= SOL_cg[i];
  for (i=0; i<nB; i++) zbsub[i] = TEMP_cg[dofB[i]]*weight[i];
  //  std::cout << *zB_Sub << std::endl;
  //
  // coarse grid correction
  //
  for (i=0; i<nB; i++) vecphib[i] = rbsub[i];
  PhiB->Multiply(true, *vec_PhiB, *prod_PhiB);
  work_Kc->Import(*prod_PhiB, *Importer_Kc, Insert);
  nKc = work_Kc->GlobalLength();
  if (MyPID == 0) {
    for (i=0; i<num_tied_down; i++) workkc[tied_down[i]] = 0;
    if (nKc > 0) AKc->sol(1, workkc, SOL_Kc, TEMP_Kc);
    for (i=0; i<nKc; i++) workkc[i] = SOL_Kc[i];
  }
  prod_PhiB->Export(*work_Kc, *Exporter_Kc, Insert);
  PhiB->Multiply(false, *prod_PhiB, *vec_PhiB);
  for (i=0; i<nB; i++) zbsub[i] += vecphib[i]*weight[i];
  //  std::cout << *zB_Sub << std::endl;
  zB_St->Export(*zB_Sub, *ExporterB, Add);
  //  std::cout << *zB_St << std::endl;
  //  if (num_tied_down > 0) remove_orthog_null(zB_St);
}

void CLIP_solver::final_update(Epetra_Vector* uB, Epetra_Vector* rB)
{
  int i;
  double tpe, *ub, *rb;
  Epetra_BLAS EB;
  char TRANS = 'T'; double ALPHA(1); double BETA(0);
  uB->ExtractView(&ub);
  rB->ExtractView(&rb);
  if ((nB_own > 0) && (n_orthog_used > 0)) 
    EB.GEMV(TRANS, nB_own, n_orthog_used, ALPHA, P_ortho, nB_own, rb, BETA, 
	    orth1);
  else
    myzero(orth1, n_orthog_used);
  Comm.SumAll(orth1, orth2, n_orthog_used);
  TRANS = 'N'; BETA = 1;
  if ((nB_own > 0) && (n_orthog_used > 0)) {
    EB.GEMV(TRANS, nB_own, n_orthog_used, ALPHA, P_ortho, nB_own, orth2, BETA, 
	    ub);
  }
}

void CLIP_solver::remove_search(Epetra_Vector* v)
{
  int j, M, N, LDA, number_cgs(2);
  char TRANS;
  double ALPHA, BETA, *AP, *P, *X, *Y;
  Epetra_BLAS EB;
  AP = AP_ortho;
  P  =  P_ortho;
  v->ExtractView(&X); 
  M = zB_St->MyLength();
  N = n_orthog_used;
  LDA = M;
  for (j=0; j<number_cgs; j++) {
    TRANS = 'T'; ALPHA = 1; BETA = 0;
    if ((M > 0) && (N > 0)) 
      EB.GEMV(TRANS, M, N, ALPHA, AP, LDA, X, BETA, orth1);
    else 
      myzero(orth1, N);
    Comm.SumAll(orth1, orth2, N);
    TRANS = 'N'; ALPHA = -1; BETA = 1;
    if ((M > 0) && (N > 0)) 
      EB.GEMV(TRANS, M, N, ALPHA,  P, LDA, orth2, BETA, X);
  }
}

void CLIP_solver::determine_AxB(Epetra_Vector* x, Epetra_Vector* b)
{
  int i, j, NumEntries, *Indices, dof;
  double *Values, sum, *pbsub, *rbsub, dnorm;
  pB_Sub->Import(*x, *ImporterB, Insert);
  pB_Sub->ExtractView(&pbsub);
  rB_Sub->ExtractView(&rbsub);
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nB; i++) {
    dof = dofB[i];
    Block_Stiff->ExtractMyRowView(dof, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) TEMP_cg[Indices[j]] -= Values[j]*pbsub[i];
  }
  for (i=0; i<nI; i++) RHS_cg[i] = TEMP_cg[dofI[i]];
  if (nI > 0) AI->sol(1, RHS_cg, SOL_cg, TEMP_cg);
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nB; i++) TEMP_cg[dofB[i]] = pbsub[i];
  for (i=0; i<nI; i++) TEMP_cg[dofI[i]] = SOL_cg[i];
  myzero(rbsub, nB);
  for (i=0; i<nB; i++) {
    dof = dofB[i];
    Block_Stiff->ExtractMyRowView(dof, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) rbsub[i] += Values[j]*TEMP_cg[Indices[j]];
  }
  b->PutScalar(0.0);
  b->Export(*rB_Sub, *ExporterB, Add);
  /*
  work_Sub->PutScalar(0.0);
  double *worksub;
  work_Sub->ExtractView(&worksub);
  for (i=0; i<ndof_sub; i++) {
    Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
    sum = 0;
    for (j=0; j<NumEntries; j++) sum += Values[j]*TEMP_cg[Indices[j]];
    worksub[i] = sum;
  }
  work_St->Export(*work_Sub, *Exporter, Add);
  b->Norm2(&dnorm);
  if (MyPID == 0) std::cout << "2 norm of ApB_St = " << dnorm << std::endl;
  work_St->Norm2(&dnorm);
  if (MyPID == 0) std::cout << "2 norm of Ap_St = " << dnorm << std::endl;
  */
}

void CLIP_solver::store_search(double pAp)
{
  int i, ibeg;
  double sfac, *pbst, *apbst;
  if (n_orthog_used < max_orthog) {
    pB_St->ExtractView(&pbst);
    ApB_St->ExtractView(&apbst);
    sfac = 1/sqrt(pAp);
    ibeg = n_orthog_used*nB_own;
    for (i=0; i<nB_own; i++) {
      P_ortho[ibeg]  = sfac*pbst[i];
      AP_ortho[ibeg] = sfac*apbst[i];
      ibeg++;
    }
    n_orthog_used++;
  }
  if ((n_orthog_used == max_orthog) && (cg_iter == -1)) cg_iter = 0;
}

void CLIP_solver::calculate_u_St()
{
  int i, j, NumEntries, *Indices, dof;
  double *Values, sum, *ubsub, *rsub, *usub;
  Epetra_Vector *uB_Sub = rB_Sub;
  Epetra_Vector *u_Sub = work_Sub;

  uB_Sub->Import(*uB_St, *ImporterB, Insert);
  uB_Sub->ExtractView(&ubsub);
  r_Sub->Import(*r_St, *Importer, Insert);
  r_Sub->ExtractView(&rsub);
  u_Sub->ExtractView(&usub);

  for (i=0; i<nB; i++) {
    dof = dofB[i];
    Block_Stiff->ExtractMyRowView(dof, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) rsub[Indices[j]] -= Values[j]*ubsub[i];
  }
  for (i=0; i<nI; i++) RHS_cg[i] = rsub[dofI[i]];
  if (nI > 0) AI->sol(1, RHS_cg, SOL_cg, TEMP_cg);
  myzero(usub, ndof_sub);
  for (i=0; i<nI; i++) usub[dofI[i]] = SOL_cg[i];
  for (i=0; i<nB; i++) usub[dofB[i]] = ubsub[i];
  u_St->Export(*u_Sub, *Exporter, Insert);
}


void CLIP_solver::calculate_Au(Epetra_Vector *x_St, Epetra_Vector *Ax_St)
{
  int i, j, NumEntries, *Indices, dof;
  double *Values, sum, *xsub, *axsub;
  Epetra_Vector *x_Sub = work_Sub;
  Epetra_Vector *Ax_Sub = r_Sub;

  x_Sub->Import(*x_St, *Importer, Insert);
  x_Sub->ExtractView(&xsub);
  Ax_Sub->ExtractView(&axsub);
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<ndof_sub; i++) {
    Block_Stiff->ExtractMyRowView(i, NumEntries, Values, Indices);
    sum = 0;
    for (j=0; j<NumEntries; j++) sum += Values[j]*xsub[Indices[j]];
    axsub[i] = sum;
  }
  Ax_St->Export(*Ax_Sub, *Exporter, Add);
}

void CLIP_solver::calculate_AuStand(Epetra_Vector *u, Epetra_Vector *Au)
{
  int i, j, NumEntries, *Indices, nStand;
  double *Values, sum, *vec1, *vec2;

  nStand = ASub->NumMyRows();
  vec1_ASub->Import(*u, *ImporterStand, Insert);
  vec1_ASub->ExtractView(&vec1); 
  vec2_ASub->ExtractView(&vec2); myzero(vec2, nStand);
  for (i=0; i<nStand; i++) {
    ASub->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) vec2[i] += Values[j]*vec1[Indices[j]];
  }
  Au->Export(*vec2_ASub, *ExporterStand, Add);
}

void CLIP_solver::two_steps_CGS(int gmres_iter, Epetra_Vector* r)
{
  int i;
  double dprod, normr;
  Epetra_BLAS EB;
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
      std::cout << "orthogonality for i = " << i << std::endl;
      for (int j=0; j<N+1; j++) std::cout << gmres_sum[j] << std::endl;
    }
  }
  */
}

void CLIP_solver::hessenberg_qr(int gmres_iter)
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

void CLIP_solver::gmres_givens(double a, double b, double & c, double & s)
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

void CLIP_solver::construct_solution(int gmres_iter, double normb)
{
  int i, j, ii, n;
  double w1, w2;
  Epetra_BLAS EB;
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
  double *Y = rB_St->Values();
  int M = rB_St->MyLength();
  int N = gmres_iter;
  int LDA = M;
  double ALPHA(normb); double BETA(0);
  if (M > 0) EB.GEMV(TRANS, M, N, ALPHA, A, LDA, X, BETA, Y);
  apply_preconditioner();
  uB_St->Update(1.0, *zB_St, 0.0);
}

void CLIP_solver::calculate_multipliers(Epetra_Vector* uStand, 
                double & norm_rconstraint, double & norm_conerror)
{
  int i;
  double *Values;
  vStand->ExtractView(&Values);
  for (i=0; i<nmycon; i++) {
    lambda_local[i] = 0;
    if (mycdof[i] >= 0) lambda_local[i] = Values[mycdof[i]];
  }
  Lambda->Export(*Lambda_local, *Exporter_lam, Insert);
  CtT->Multiply(false, *Lambda, *wStand);
  vStand->Update(-1.0, *wStand, 1.0); 
  vStand->Norm2(&norm_rconstraint);
  ConStandard->Multiply(true, *uStand, *ConError);
  double norm_uStand, norm_ConError, norm_ConStandard;
  uStand->NormInf(&norm_uStand);
  ConError->NormInf(&norm_ConError);
  norm_ConStandard = ConStandard->NormInf();
  norm_conerror = norm_ConError/norm_ConStandard/norm_uStand;
  //  Lambda->Export(*Lambda_local, *Exporter_lam, Insert);  
}

void CLIP_solver::mpcforces(Epetra_Vector* uLocal, 
			    Epetra_Import* ImporterStLam2Loc)
{
}

void CLIP_solver::calculate_condition(int miter)
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
    DSTEV_F77(&N, &ip1, Dtri, Etri, &Z, &one, &WORK, &INFO, 1); 
    if (INFO != 0) {
      std::cout << "error in call to DSTEV in CLASP_solver::calculate_condition" 
	   << std::endl;
      std::cout << "INFO = " << INFO << std::endl;
    }
    econa[i] = Dtri[i]/Dtri[0];
  }

}

void CLIP_solver::spmat_datfile(const Epetra_CrsMatrix & A, char fname[], 
				int opt)
{
  int i, j, NumEntries, *Indices, grow, gcol;
  double *Values;
  std::ofstream ffout;
  ffout.open(fname);
  for (i=0; i<A.NumMyRows(); i++) { 
    A.ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (opt == 1)
	ffout << i+1 << " " << Indices[j]+1 << std::setw(22) << std::setprecision(15)
	     << Values[j] << std::endl;
      if (opt == 2) {
	grow = A.GRID(i); gcol = A.GCID(Indices[j]);
	ffout << grow+1 << " " << gcol+1 << std::setw(22) 
	     << std::setprecision(15) << Values[j] << std::endl;
      }
    }
  }
  ffout.close();
}
