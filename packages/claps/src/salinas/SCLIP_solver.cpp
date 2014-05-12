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
#include <assert.h>
#include <algorithm>
#include "SCLIP_solver.hpp"
#include "myzero.hpp"

SCLIP_solver::SCLIP_solver(
  int nnode_,               // number of nodes in superstructure
  int nelem_,               // number of elements in superstructure
  int neq_,                 // number of of active dofs in superstructure
  int nadj_proc_,           // number of adjacent processors
  int ndofs_deleted_,       // number of deleted degrees of freedom
  int NumberMpc_,           // number of constraint equations
  const int E1_[],          // nodes (local numbering) for elements
  const int E2_[],          // pointer array for E1_
  const int adj_proc_[],    // adjacent processors
  const int H1_global_[],   // nodes on boundary (global numbering)
  const int H2_[],          // pointer array for H1_global
  const MPI_Comm *mpicomm_, // pointer to MPI communicator
  const int gnn_[],         // global node numbers in superstructure
  const unsigned short *dofmap_on_node_, // codes for active nodal dofs
  const double x_[],        // x-coordinates of nodes
  const double y_[],        // y-coordinates of nodes
  const double z_[],        // z-coordinates of nodes
  const int rowbeg_KU_[],   // beginning of rows for      upper t sparse matrix
  const int colidx_KU_[],   // nonzero column numbers in  upper t sparse matrix
  const double KU_[],       // nonzero entries        in  upper t sparse matrix
  const int mapdof2row_[],  // a degree of freedom map
  const int bc_node_ids_[], // array of deleted dof node numbers
  const int bc_node_dof_[], // array of deleted dof local dof numbers
  const MpcLocal* MpcVector_)
{
  nnode          = nnode_;
  nelem          = nelem_;
  neq            = neq_;
  nadj_proc      = nadj_proc_;
  ndofs_deleted  = ndofs_deleted_;
  NumberMpc      = NumberMpc_;
  E1             = E1_;
  E2             = E2_;
  adj_proc       = adj_proc_;
  H1_global      = H1_global_;
  H2             = H2_;
  mpicomm        = *mpicomm_;
  gnn            = gnn_;
  dofmap_on_node = dofmap_on_node_;
  x              = x_;
  y              = y_;
  z              = z_;
  rowbeg_KU      = rowbeg_KU_;
  colidx_KU      = colidx_KU_;
  KU             = KU_;
  mapdof2row     = mapdof2row_;
  bc_node_ids    = bc_node_ids_;
  bc_node_dof    = bc_node_dof_;
  MpcVector      = MpcVector_;
  Comm = new Epetra_MpiComm(mpicomm);
  MyPID = Comm->MyPID();
  maxdofnode     = 7;
  ASub = 0; ConStandard = 0; Coords = 0; SubMap = 0; StandardMap = 0;
  NodalDofs = 0; H1 = 0; uStand = 0; uLocal = 0;
  fStand = 0; uSub = 0; uvec = 0; fvec = 0; subvec = 0; locvec = 0;
  ImporterSt2Sub = 0; ExporterSub2St = 0; MpcLocalMap = 0;
  ImporterStLam2Loc = 0; CS = 0;
}

SCLIP_solver::~SCLIP_solver()
{
  delete Comm; delete ASub; delete ConStandard; delete Coords;
  delete NodalDofs; delete uStand; delete uLocal;
  delete fStand; delete uSub; delete [] uvec; delete [] fvec;
  delete [] subvec; delete [] locvec; delete SubMap; delete StandardMap; 
  delete ImporterSt2Sub; delete ExporterSub2St; delete MpcLocalMap;
  delete ImporterStLam2Loc; delete CS;
}

void SCLIP_solver::construct_H1()
{
  int i, j, n;
  H1 = 0;
  if (nnode > 0) {
    //
    // assign sorted H1_global segments to H1
    //
    H1 = new int[H2[nadj_proc]];
    for (i=0; i<nadj_proc; i++) {
      for (j=H2[i]; j<H2[i+1]; j++) H1[j] = H1_global[j];
      n = H2[i+1] - H2[i];
      if (n > 1) std::sort(&H1[H2[i]], &H1[H2[i+1]]);
    }
    //
    // sort global node numbers
    //
    int *gnn_sort = new int[nnode];
    for (i=0; i<nnode; i++) gnn_sort[i] = gnn[i];
    std::sort(&gnn_sort[0], &gnn_sort[nnode]);
    int *sort_to_orig = new int[nnode];
    for (i=0; i<nnode; i++) {
      int sorted_loc = CRD_utils::find_index(gnn_sort, nnode, gnn[i]);
      assert(sorted_loc != -1);
      sort_to_orig[sorted_loc] = i;
    }
    //
    // convert H1 global node numbers to local node numbers
    //
    for (i=0; i<H2[nadj_proc]; i++) {
      int sorted_loc = CRD_utils::find_index(gnn_sort, nnode, H1[i]);
      assert(sorted_loc != -1);
      H1[i] = sort_to_orig[sorted_loc];
    }
    delete [] gnn_sort; delete [] sort_to_orig;
  }
}

void SCLIP_solver::determine_ownership()
{
  int i, j, k, m, n, node, isum, dof, sproc, rproc, ierr;
  int a[7] = {1, 2, 4, 8, 16, 32, 64};
  //
  // determine active degrees of freedom
  //   nodedof[nodebeg[i]:nodebeg[i+1]-1] = active dofs for node i
  //   sub_gdof[0:ndof-1] = active global dofs for processor
  //
  unsigned char *active_dof = new unsigned char[maxdofnode*nnode];
  myzero(active_dof, maxdofnode*nnode);
  for (i=0; i<ndofs_deleted; i++) {
    active_dof[maxdofnode*bc_node_ids[i]+bc_node_dof[i]] = 1;
  }
  ndof = 0;
  for (i=0; i<nnode; i++) {
    m = maxdofnode*i;
    for (j=0; j<maxdofnode; j++) {
      if ((dofmap_on_node[i] &  a[j]) && (active_dof[m+j] == 0)) ndof++;
    }
  }
  int *sub_gdof = new int[ndof];
  int *nodebeg = new int[nnode+1]; nodebeg[0] = 0;
  unsigned char *nodedof = new unsigned char[ndof];
  ndof = 0;
  for (i=0; i<nnode; i++) {
    m = maxdofnode*i;
    for (j=0; j<maxdofnode; j++) {
      if ((dofmap_on_node[i] &  a[j]) && (active_dof[m+j] == 0)) {
	nodedof[ndof] = j;
	sub_gdof[ndof] = maxdofnode*gnn[i] + j;
	ndof++;
      }
    }
    nodebeg[i+1] = ndof;
  }
  assert(ndof == neq);
  SubMap = new Epetra_Map(-1, ndof, sub_gdof, 0, *Comm);
  delete [] active_dof;
  //
  // constraints
  //
  //
  // first determine number of null constraint equations
  //
  int ncon_zero(0);
  for (i=0; i<NumberMpc; i++) 
    if (MpcVector[i].NumEntriesGlobal() == 0) ncon_zero++;
  if ((MyPID == 0) && (ncon_zero > 0)) {
    std::cout << "Warning: number of null constraints = " << ncon_zero;
  }
  int *count1 = new int[ndof]; myzero(count1, ndof);
  int *flag_con = new int[NumberMpc]; myzero(flag_con, NumberMpc);
  int *flag_con_max = new int[NumberMpc];
  nmpc_loc = 0;
  for (i=0; i<NumberMpc; i++) {
    if (MpcVector[i].NumEntries() > 0) nmpc_loc++;
    for (j=0; j<MpcVector[i].NumEntries(); j++) {
      node = MpcVector[i].LocalId(j);
      dof  = MpcVector[i].NodalDof(j);
      for (k=nodebeg[node]; k<nodebeg[node+1]; k++) {
	if (nodedof[k] == dof) {
	  flag_con[i] = 1;
	  count1[k]++;
	  break;
	}
      }
    }
  }
  Comm->MaxAll(flag_con, flag_con_max, NumberMpc);
  delete [] flag_con;
  //
  // determine number of constraints with no active dofs (each of these 
  // constraints will be ignored)
  //
  if (MyPID == 0) {
    int ncon_no_active(0);
    for (i=0; i<NumberMpc; i++) if (flag_con_max[i] == 0) ncon_no_active++;
    if (ncon_no_active > 0) {
      std::cout << "Warning: number of constraints with no active dofs = " 
	        << ncon_no_active << std::endl;
      std::cout << "These constraint equations will be ignored" << std::endl;
    }
  }
  int *con2 = new int[ndof+1]; con2[0] = 0;
  for (i=0; i<ndof; i++) con2[i+1] = con2[i] + count1[i];
  myzero(count1, ndof);
  int *con1 = new int[con2[ndof]];
  double *coef = new double[con2[ndof]]; myzero(coef, con2[ndof]);
  int *mpc_loc = new int[nmpc_loc]; 
  locvec = new double[nmpc_loc]; nmpc_loc = 0;
  int NumberMpc_keep(0);
  for (i=0; i<NumberMpc; i++) {
    if (flag_con_max[i] > 0) {
      for (j=0; j<MpcVector[i].NumEntries(); j++) {
	node = MpcVector[i].LocalId(j);
	dof  = MpcVector[i].NodalDof(j);
	for (k=nodebeg[node]; k<nodebeg[node+1]; k++) {
	  if (nodedof[k] == dof) {
	    con1[con2[k] + count1[k]] = NumberMpc_keep;
	    coef[con2[k] + count1[k]] = MpcVector[i].Coef(j);
	    count1[k]++;
	    break;
	  }
	}
      }
      if (MpcVector[i].NumEntries() > 0) {
	mpc_loc[nmpc_loc] = NumberMpc_keep;
	nmpc_loc++;
      }
      if (MpcVector[i].NumEntriesGlobal() > 0) NumberMpc_keep++;
    }
  }
  delete [] flag_con_max;
  Epetra_Map RowMapCon(NumberMpc_keep, 0, *Comm);
  MpcLocalMap = new Epetra_Map(-1, nmpc_loc, mpc_loc, 0, *Comm);
  Epetra_CrsMatrix ConSubdomain(View, *SubMap, count1);
  for (i=0; i<nnode; i++) {
    for (j=nodebeg[i]; j<nodebeg[i+1]; j++) {
      int gdof = maxdofnode*gnn[i] + nodedof[j];
      ierr = ConSubdomain.InsertGlobalValues(gdof, count1[j], &coef[con2[j]],
					     &con1[con2[j]]);
      assert (ierr == 0);
    }
  }
  delete [] mpc_loc;
  //
  // determine dof_ownership
  //
  int *H1codeS = new int[H2[nadj_proc]];
  int *H1codeR = new int[H2[nadj_proc]];
  for (i=0; i<H2[nadj_proc]; i++) {
    node = H1[i];
    isum = 0;
    for (j=nodebeg[node]; j<nodebeg[node+1]; j++) {
      dof = nodedof[j];
      isum += a[dof];
    }
    H1codeS[i] = isum;
  }
  MPI_Request *arequest = new MPI_Request[nadj_proc];
  MPI_Status *astatus  = new MPI_Status[nadj_proc];
  int tag = 7;
  for (i=0; i<nadj_proc; i++) {
    sproc = adj_proc[i];
    n = H2[i+1] - H2[i];
    ierr = MPI_Irecv(&H1codeR[H2[i]], n, MPI_INT, sproc, tag, 
		     mpicomm, &arequest[i]);
  }
  for (i=0; i<nadj_proc; i++) {
    rproc = adj_proc[i];
    n = H2[i+1] - H2[i];
    ierr = MPI_Send( &H1codeS[H2[i]], n, MPI_INT, rproc, tag, mpicomm);
  }
  ierr = MPI_Waitall(nadj_proc, arequest, astatus);
  assert (ierr == 0);
  delete [] H1codeS; delete [] arequest; delete [] astatus;
  bool *iown = new bool[maxdofnode*nnode];
  for (i=0; i<maxdofnode*nnode; i++) iown[i] = true;
  for (i=0; i<nadj_proc; i++) {
    for (j=H2[i]; j<H2[i+1]; j++) {
      m = maxdofnode*H1[j];
      for (k=0; k<maxdofnode; k++) {
	if ((H1codeR[j] &  a[k]) && (adj_proc[i] < MyPID)) iown[m+k] = false;
      }
    }
  }
  delete [] H1codeR; delete [] H1;
  ndof_mine = 0;
  for (i=0; i<nnode; i++) {
    m = maxdofnode*i;
    for (j=nodebeg[i]; j<nodebeg[i+1]; j++) {
      dof = nodedof[j];
      if (iown[m+dof] == true) {
	sub_gdof[ndof_mine] = maxdofnode*gnn[i] + dof;
	ndof_mine++;
      }
    }
  }
  uvec = new double[ndof_mine]; fvec = new double[ndof_mine];
  //
  // coordinates, local dofs, and local nodes for substructure dofs
  //
  double *xyzdof = new double[3*ndof];
  int *localdofs = new int[ndof];
  int ndof_mine_keep = ndof_mine;
  int jlow = 0;
  for (i=0; i<nnode; i++) {
    for (j=nodebeg[i]; j<nodebeg[i+1]; j++) {
      dof = nodedof[j];
      xyzdof[j         ] = x[i];
      xyzdof[j +   ndof] = y[i];
      xyzdof[j + 2*ndof] = z[i];
      localdofs[j]       = dof + 1;
    }
  }
  delete [] iown; delete [] nodebeg; delete [] nodedof; 

  StandardMap = new Epetra_Map(-1, ndof_mine, sub_gdof, 0, *Comm);
  ConStandard = new Epetra_CrsMatrix(Copy, *StandardMap, 0);
  ierr = ConSubdomain.FillComplete(RowMapCon, *StandardMap);  
  assert (ierr == 0);
  Epetra_Export Exporter(*SubMap, *StandardMap);
  ierr = ConStandard->Export(ConSubdomain, Exporter, Add);
  assert (ierr == 0);
  //  ConStandard->Export(ConSubdomain, Exporter, Insert);
  ierr = ConStandard->FillComplete(RowMapCon, *StandardMap);
  assert (ierr == 0);
  CRD_utils::scale_columns(ConStandard, 1, 1000);
  //  cout << *ConStandard << endl;

  delete [] count1; delete [] con1; delete [] con2; delete [] coef;
  Coords = new Epetra_MultiVector(Copy, *SubMap, xyzdof, ndof, 3);
  NodalDofs = new Epetra_IntVector(Copy, *SubMap, localdofs);
  uStand = new Epetra_Vector(View, *StandardMap, uvec);
  fStand = new Epetra_Vector(View, *StandardMap, fvec);
  uLocal = new Epetra_Vector(View, *MpcLocalMap, locvec);
  delete [] sub_gdof; delete [] xyzdof; delete [] localdofs; 
}

void SCLIP_solver::construct_K_base()
{
  int i, j, dofrow, dofcol, ierr, ibeg, ndiag(0);
  int *count = new int[neq];
  myzero(count, neq);
  int *maprow2dof = new int[neq];
  for (i=0; i<neq; i++) maprow2dof[mapdof2row[i]] = i;
  for (i=0; i<neq; i++) {
    dofrow = maprow2dof[i];
    for (j=(rowbeg_KU[i]-1); j<(rowbeg_KU[i+1]-1); j++) {
      if ((colidx_KU[j]-1) != i) {
	dofcol = maprow2dof[(colidx_KU[j]-1)];
	count[dofcol]++;
      }
      else ndiag++;
      count[dofrow]++;
    }
  }
  int *rowbeg_base = new int[neq+1]; rowbeg_base[0] = 0;
  for (i=0; i<neq; i++) rowbeg_base[i+1] = rowbeg_base[i] + count[i];
  int nnz = rowbeg_base[neq];
  int nnz_KU =rowbeg_KU[neq]-1;
  if (neq == 0) nnz_KU = 0;
  assert(nnz == (2*nnz_KU-ndiag));
  int *colidx_base = new int[nnz];
  double *K_base = new double[nnz];
  myzero(count, neq);
  for (i=0; i<neq; i++) {
    dofrow = maprow2dof[i];
    for (j=(rowbeg_KU[i]-1); j<(rowbeg_KU[i+1]-1); j++) {
      dofcol = maprow2dof[colidx_KU[j]-1];
      if ((colidx_KU[j]-1) != i) {
	colidx_base[rowbeg_base[dofcol] + count[dofcol]] = dofrow;
	K_base[     rowbeg_base[dofcol] + count[dofcol]] = KU[j];
	count[dofcol]++;
      }
      colidx_base[rowbeg_base[dofrow] + count[dofrow]] = dofcol;
      K_base     [rowbeg_base[dofrow] + count[dofrow]] = KU[j];
      count[dofrow]++;
    }
  }
  delete [] maprow2dof;
  ASub = new Epetra_CrsMatrix(Copy, *SubMap, *SubMap, count);
  for (i=0; i<neq; i++) {
    ibeg = rowbeg_base[i];
    ierr = ASub->InsertMyValues(i, count[i], &K_base[ibeg], 
				&colidx_base[ibeg]);
  }
  ASub->FillComplete();
  ASub->OptimizeStorage();

  Epetra_Export Exporter(*SubMap, *StandardMap);

  int NumEntries, *Indices;
  double min_diag_all, max_diag_all, *Values;
  Epetra_Vector ASub_diag(*SubMap);
  Epetra_Vector ASt_diag(*StandardMap);
  ASub_diag.PutScalar(0.0);
  for (i=0; i<ASub->NumMyRows(); i++) {
    ASub->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (Indices[j] == i) ASub_diag[i] = Values[j];
    }
  }
  ASt_diag.Export(ASub_diag, Exporter, Add);
  ASt_diag.MinValue(&min_diag_all);
  ASt_diag.MaxValue(&max_diag_all);
  if (MyPID == 0) {
    //    cout << "min diagonal entry in AStandard = " << min_diag_all << endl;
    //    cout << "max diagonal entry in AStandard = " << max_diag_all << endl;
    if (min_diag_all != 0) { 
      //      cout << "ratio of max to min diagonal = "  
      //	   << max_diag_all/min_diag_all << endl;
    }
  }
  int print_flag = 1;
  if (print_flag == 0) {
    char filename[101];
    sprintf(filename,"%s%d","CLIP_base", MyPID);
    CRD_utils::spmat_datfile(neq ,rowbeg_base, colidx_base, K_base, filename);
    EPmat_datfile(ASub, filename);
  }

  delete [] rowbeg_base; delete [] colidx_base; delete [] K_base; 
  delete [] count;
  ImporterSt2Sub = new Epetra_Import(*SubMap, *StandardMap);
  ExporterSub2St = new Epetra_Export(*SubMap, *StandardMap);
  ImporterStLam2Loc = new Epetra_Import(*MpcLocalMap, 
					ConStandard->DomainMap());
  subvec = new double[ndof];
  uSub = new Epetra_Vector(View, *SubMap, subvec);
}

void SCLIP_solver::CLIP_solver_init(int cdof_option, double solver_tol,
   int maxiter, int atype, int ndim, int local_solver, int max_orthog,
   int prt_debug, int prt_summary, int chk_sub_singularity, 
   int krylov_method, int scale_option, int num_rigid_mode)
{
  double clip_params[20];

  clip_params[0]  = double(cdof_option);
  clip_params[1]  = solver_tol;
  clip_params[2]  = double(maxiter);
  clip_params[3]  = double(max_orthog);
  clip_params[4]  = double(atype);
  clip_params[5]  = double(ndim);
  clip_params[6]  = double(local_solver);
  clip_params[7]  = double(prt_debug);
  clip_params[8]  = double(prt_summary);
  clip_params[9]  = double(chk_sub_singularity);
  clip_params[10] = double(krylov_method);
  clip_params[11] = double(scale_option);
  clip_params[12] = double(num_rigid_mode);
  CS = new CLIP_solver(ASub, NodalDofs, Coords, ConStandard, clip_params);
}

void SCLIP_solver::solve(double f[], double u[], int & number_iterations, 
			 int & SCLIP_status, int & max_added_corner)
{
  int i;
  for (i=0; i<ndof; i++) subvec[i] = f[i];
  fStand->Export(*uSub, *ExporterSub2St, Add);
  CS->solve(uStand, fStand, number_iterations, SCLIP_status, max_added_corner);
  uSub->Import(*uStand, *ImporterSt2Sub, Insert);
  for (i=0; i<ndof; i++) u[i] = subvec[i];
}

void SCLIP_solver::MpcForces( double *cvals)
{
  int i;
  //  CS->mpcforces(uLocal, ImporterStLam2Loc);
  //  for (i=0; i<nmpc_loc; i++) cvals[i] = locvec[i];
  for (i=0; i<nmpc_loc; i++) cvals[i] = 0;
}

void SCLIP_solver::EPmat_datfile(Epetra_CrsMatrix* A, char fname[])
{
  int i, j, NumEntries, *Indices;
  double *Values;
  std::ofstream fout;
  sprintf(fname, "%s.epetra", fname);
  fout.open(fname);
  for (i=0; i<A->NumMyRows(); i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      fout << A->GRID(i)+1 << " " << A->GCID(Indices[j])+1 << std::setw(22) 
	   << std::setprecision(15) << Values[j] << std::endl;
    }
  }
  fout.close();
}
