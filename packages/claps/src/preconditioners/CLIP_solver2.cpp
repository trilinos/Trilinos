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

#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <assert.h>
//#include <climits>
#include "CLIP_solver2.hpp"
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

CLIP_solver2::CLIP_solver2(CRS_serial* A_,
			   const Epetra_Map* SubMap_,
			   const Epetra_Map* OwnMap_,
			   const double* clip_params_,
			   const double* amg_params_)
  : A(A_), SubMap(SubMap_), OwnMap(OwnMap_), clip_params(clip_params_),
    amg_params(amg_params_), Comm(SubMap->Comm())
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
  sub_solver          = int(clip_params[13]);
  coarse_solver       = int(clip_params[14]);

  double starttime, endtime;
  Comm.Barrier();
  MyPID = Comm.MyPID();
  if (MyPID == 0) starttime = MPI_Wtime();
  NumProc = Comm.NumProc();
  ndof_global = OwnMap->NumGlobalElements();
  print_flag = -1;
  if (MyPID == 0) {
    print_flag = prt_summary + 10*prt_debug;
    fout.open("CLIP_solver2.data");
    fout << "----------------- CLIP solver summary information "
	 << "-----------------" << endl;
  }
  if (print_flag > 0) fout << "number of global dofs        = " << ndof_global 
			   << endl;
  Exporter = new Epetra_Export(*SubMap, *OwnMap);
  Importer = new Epetra_Import(*SubMap, *OwnMap);
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
  // allocate memory for vectors used in preconditioner
  //
  allocate_vectors();
  //
  //  const Epetra_MpiComm &empicomm = 
  //               dynamic_cast<const Epetra_MpiComm &>(Comm);
  //  mpicomm = empicomm.Comm();
  Comm.Barrier();
  if (print_flag > 0) {
    endtime = MPI_Wtime();
    fout << "elapsed time for CLIP solver init  = " 
	 << endtime-starttime << " seconds" << endl;
    if (krylov_method != 1) fout << "pcg solver will be used" << endl;
    if (krylov_method == 1) fout << "gmres solver will be used" << endl;
    fout.close();
  }
}

CLIP_solver2::~CLIP_solver2()
{
  delete Importer; delete ImporterB; delete ImporterC;
  delete Exporter; delete ExporterB; delete ExporterC;
  delete SubMapB; delete SubMapC; delete OwnMapB; delete OwnMapC;
  delete [] PhiB; 
  delete AR; delete AI; delete AKc;
  delete [] comp1; delete [] comp2; delete [] sub1; delete [] sub2;
  delete [] dset1; delete [] dset2; delete [] corner_flag; 
  delete [] dofI; delete [] dofB; delete [] dofC; delete [] dofR;
  delete [] dofBB; delete [] weight; delete [] dof2node;
  delete [] ARinvCT; delete [] CARinvCT; delete [] lambda_r; delete [] RHS_cg; 
  delete [] SOL_cg; delete [] TEMP_cg; delete [] SOL_Kc; delete [] TEMP_Kc;
  delete [] amg_matR; delete [] amg_A1R; delete [] amg_A2R; delete [] amg_xR;
  delete [] amg_yR; delete [] amg_zR; delete [] amg_nodebegR; delete [] amg_local_dofR;
  delete [] amg_matI; delete [] amg_A1I; delete [] amg_A2I; delete [] amg_xI;
  delete [] amg_yI; delete [] amg_zI; delete [] amg_nodebegI; delete [] amg_local_dofI;
  delete [] amg_matC; delete [] amg_A1C; delete [] amg_A2C; delete [] amg_xC;
  delete [] amg_yC; delete [] amg_zC; delete [] amg_nodebegC; delete [] amg_local_dofC;
  delete rSub; delete uSub ; delete rOwn; delete uOwn; delete rSubB; delete uSubB;
  delete rOwnB; delete uOwnB; delete rSubC; delete uSubC; delete rOwnC; delete uOwnC;
}

void CLIP_solver2::determine_components()
{
  //
  // determine subdomain components:
  //   ncomp = number of components
  //   comp1[comp2[i]:comp2[i+1]-1] = subdomain dofs in component i
  //
  int *a1 = A->get_colidx();
  int *a2 = A->get_rowbeg();
  determine_components(a1, a2, ndof_sub, comp1, comp2, ncomp);
  //  cout << "MyPID, ncomp = " << MyPID << " " << ncomp << endl;
}

void CLIP_solver2::determine_dof_sets()
{
  int i, j, k, NumEntries, *Indices, dof, dof2, *loc_dofs;
  int ScanSums;
  loc_dofs = A->get_localdof();
  //
  // dof2node[i] = local node number for dof i
  //
  dof2node = new int[ndof_sub];
  int* nodebeg = A->get_nodebeg();
  int nnode = A->get_nnode();
  for (i=0; i<nnode; i++)
    for (j=nodebeg[i]; j<nodebeg[i+1]; j++) dof2node[j] = i;
  //
  // each subdomain component is now treated and numbered as an
  // individual subdomain
  //  ncomp_sum = total number of subdomains
  //
  Comm.ScanSum(&ncomp, &ScanSums, 1);
  Comm.SumAll(&ncomp, &ncomp_sum, 1);
  //  cout << "MyPID, ScanSums = " << MyPID << " " << ScanSums << endl;
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
  Subdomains.FillComplete(ColMap, *OwnMap);
  Epetra_CrsGraph SubStand(Copy, *OwnMap, 0);
  SubStand.Export(Subdomains, *Exporter, Insert);
  SubStand.FillComplete(ColMap, *OwnMap);
  Epetra_CrsGraph SubsForDofs(Copy, *SubMap, 0);
  SubsForDofs.Import(SubStand, *Importer, Insert);
  SubsForDofs.FillComplete(ColMap, *OwnMap);
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
  //  cout << "MyPID, ndof_set = " << MyPID << " " << ndof_set << endl;
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
    double *lower = lower_bound(&values[0], &values[ndof_set], dval_orig[i]);
    ival[i] = (int) (lower - &values[0]);
    icount[ival[i]]++;
    /*
    if (MyPID == 0) cout << "val, ival = " << dval_orig[i] << " " 
    			 << ival[i] << endl;
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
  if (MyPID == -1) {
    double* x = A->get_xcoord();
    double* y = A->get_ycoord();
    for (i=0; i<ndof_set; i++) {
      cout << "dset " << i << endl;
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	cout << x[dof2node[dset1[j]]] << " " << y[dof2node[dset1[j]]] << "; "
	     << loc_dofs[dset1[j]] << " " << dset1[j] << endl;
      }
    }
  }
  //
  // check
  //
  for (i=0; i<ndof_set; i++) {
    dof = dset1[dset2[i]];
    for (j=dset2[i]+1; j<dset2[i+1]; j++) {
      dof2 = dset1[j];
      assert(loc_dofs[dof] == loc_dofs[dof2]);
      assert((sub2[dof2+1]-sub2[dof2]) == (sub2[dof+1]-sub2[dof]));
      for (k=0; k<(sub2[dof+1]-sub2[dof]); k++) {
	assert(sub1[sub2[dof]+k] == sub1[sub2[dof2]+k]);
      }
    }
  }
  //
  // sub1 -> local column indices
  //
  for (i=0; i<ndof_sub; i++) {
    SubsForDofs.ExtractMyRowView(i, NumEntries, Indices);
    for (j=0; j<NumEntries; j++) sub1[sub2[i]+j] = Indices[j];
  }
  local_sub = new int[ncomp];
  sub_flag = new bool[SubsForDofs.NumMyCols()];
  for (i=0; i<SubsForDofs.NumMyCols(); i++) sub_flag[i] = false;
  for (i=0; i<ncomp; i++) {
    local_sub[i] = SubsForDofs.LCID(sub_begin+i);
    assert (local_sub[i] >= 0);
  }
  //
  // determine internal and boundary dofs
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
  //
  // next determine ownership for each subdomain dof
  //  owner_flag = true if subdomain owns subdomain dof i, false otherwise
  //  
  owner_flag = new bool[ndof_sub];
  for (i=0; i<ndof_sub; i++) owner_flag[i] = OwnMap->MyGID(SubMap->GID(i));  
}

void CLIP_solver2::determine_corner_dofs()
{
  int i, j, num_entries, num_sub, dof, rot_flag, *gdofs, nadof, node;
  gdofs = SubMap->MyGlobalElements();
  //
  // dofBB[1:nBB] = dofs adjacent to boundary dofs
  //  (used later for static condensation corrections)
  //
  int nnode = A->get_nnode();
  int* nodebeg = A->get_nodebeg();
  int *adof = new int[ndof_sub];
  bool *adof_flag = new bool[ndof_sub];
  for (i=0; i<ndof_sub; i++) adof_flag[i] = false;
  int *anode = new int[nnode];
  bool *anode_flag = new bool[nnode];
  for (i=0; i<nnode; i++) anode_flag[i] = false;
  get_adof(-1, dofB, nB, adof, nadof, adof_flag, 1);
  nBB = nadof - nB;
  dofBB = new int[nBB];
  memcpy(dofBB, adof+nB, nBB*sizeof(int));
  //
  // first pick singleton dof sets as corner dofs
  //
  corner_flag = new int[ndof_sub]; myzero(corner_flag, ndof_sub);
  for (i=0; i<ndof_set; i++) {
    num_entries = dset2[i+1] - dset2[i];
    dof = dset1[dset2[i]];
    if (num_entries == 1) {
      corner_flag[dof] = 1;
      node = dof2node[dof];
      for (j=nodebeg[node]; j<nodebeg[node+1]; j++) corner_flag[j] = 1;
    }
  }
  //
  // next pick corner dofs for "edges"
  //
  for (i=0; i<ncomp; i++) {
    for (j=0; j<ndof_set; j++) {
      dof = dset1[dset2[j]];
      rot_flag = determine_shell(dof2node[dof]);
      num_sub = sub2[dof+1] - sub2[dof];
      if ((find_sub(dof, local_sub[i]) == true) && (num_sub > 1)) {
	if ((num_sub > 2) || (ndim == 2) || (rot_flag == 1)) {
	  check_two_nodes(j, adof, adof_flag, anode, anode_flag);
	}
      }
    }
  }
  //
  // finally, pick corner dofs for "faces"
  //    
  for (i=0; i<ncomp; i++) {
    for (j=0; j<ndof_set; j++) {
      dof = dset1[dset2[j]];
      rot_flag = determine_shell(dof2node[dof]);
      num_sub = sub2[dof+1] - sub2[dof];
      if ((find_sub(dof, local_sub[i]) == true) && (num_sub > 1)) {
	if (!(num_sub > 2) && !(ndim == 2) && !(rot_flag == 1)) {
	  check_three_nodes(j, adof, adof_flag, anode, anode_flag);
	}
      }
    }
  }
  delete [] adof; delete [] adof_flag; delete [] sub_flag; sub_flag = 0;
  delete [] anode; delete [] anode_flag; delete [] local_sub;
  //
  // share corner dofs between subdomains
  //
  share_corner_info();
  max_added_corner = 0;
  int ncorner_dof(0);
  for (i=0; i<ndof_sub; i++) if (corner_flag[i] > 0) ncorner_dof++;
  //  cout << "MyPID, ncorner_dof = " << MyPID << " " << ncorner_dof << endl;
  if (MyPID == -1) {
    double* x = A->get_xcoord(); double* y = A->get_ycoord(); double* z = A->get_zcoord();
    for (i=0; i<ndof_sub; i++)
      if (corner_flag[i] > 0) { 
	cout << "MyPID, corner, coords = " << MyPID << " " << gdofs[i] << " " 
	     << x[dof2node[i]] << " " << y[dof2node[i]] << " " 
	     << z[dof2node[i]] << endl;
      }
  }
}

void CLIP_solver2::share_corner_info()
{
  Epetra_IntVector C_sub(View, *SubMap, corner_flag);
  Epetra_IntVector C_all(*OwnMap);
  C_all.Export(C_sub, *Exporter, Add);
  C_sub.Import(C_all, *Importer, Insert);
  for (int i=0; i<ndof_sub; i++) if (corner_flag[i] > 0) corner_flag[i] = 1;
}

void CLIP_solver2::modify_dof_sets()
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

  if (MyPID == -1) {
    cout << "dset2 = ";
    for (i=0; i<=ndof_set; i++) cout << dset2[i] << " ";
    cout << endl;
    double* x = A->get_xcoord(); double* y = A->get_ycoord();
    double* z = A->get_zcoord();
    for (i=0; i<ndof_set; i++) {
      cout << "coordinates of dofs in dof_set " << i << endl;
      for (j=dset2[i]; j<dset2[i+1]; j++) {
	dof = dset1[j];
	cout << x[dof2node[dof]] << " " << y[dof2node[dof]] << " " 
	     << z[dof2node[dof]] << endl;
      }
    }
  }

}

void CLIP_solver2::factor_sub_matrices()
{
  int i, *rowbeg, *colidx;
  double *vals;
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
  //  cout << "MyPID, atype, nC = " << MyPID << " " << atype << " " << nC << endl;
  assert ((nR + nC) == (nI + nB));
  if (sub_solver == 1) {
    gen_matrix(A, nR, dofR, rowbeg, colidx, vals);
    AR = new CLAPS_sparse_lu2();
    AR->factor(vals, rowbeg, colidx, nR, scale_option);
    delete [] rowbeg; delete [] colidx; delete [] vals;
  }
#if defined(AMG_CRD)
  if (sub_solver == 2) {
    gen_amg_data(nR, dofR, amg_matR, amg_A1R, amg_A2R, amg_xR,
		 amg_yR, amg_zR, amg_nodebegR, amg_local_dofR,
		 amg_nnodeR);
    AR = new amg_solve(amg_matR, amg_A1R, amg_A2R, amg_xR, amg_yR, amg_zR,
		       amg_nodebegR, amg_local_dofR, amg_nnodeR, amg_params);
  }
#endif
  //
  // factor AI
  //
  if (sub_solver == 1) {
    gen_matrix(A, nI, dofI, rowbeg, colidx, vals);
    AI = new CLAPS_sparse_lu2();
    AI->factor(vals, rowbeg, colidx, nI, scale_option);
    delete [] rowbeg; delete [] colidx; delete [] vals;
  }
#if defined(AMG_CRD)
  if (sub_solver == 2) {
    gen_amg_data(nI, dofI, amg_matI, amg_A1I, amg_A2I, amg_xI,
		 amg_yI, amg_zI, amg_nodebegI, amg_local_dofI,
		 amg_nnodeI);
    AI = new amg_solve(amg_matI, amg_A1I, amg_A2I, amg_xI, amg_yI, amg_zI,
		       amg_nodebegI, amg_local_dofI, amg_nnodeI, amg_params);
  }
#endif
}

void CLIP_solver2::calculate_coarse()
{
  int i, j, k, ibeg, dof, *sub_to_R, INFO, gmin, base;
  double *sub_diag, *CT, sum;

  sub_to_R = new int[ndof_sub]; for (i=0; i<ndof_sub; i++) sub_to_R[i] = -1;
  for (i=0; i<nR; i++) sub_to_R[dofR[i]] = i;
  //
  // specify global numbers for columns of Phi
  //
  int *columns = new int[nC+ndof_set];
  for (i=0; i<nC; i++) columns[i] = SubMap->GID(dofC[i]);
  base = OwnMap->MaxAllGID() + 1;
  for (i=0; i<ndof_set; i++) {
    gmin = base;
    for (j=dset2[i]; j<dset2[i+1]; j++) {
      dof = dset1[j];
      if (SubMap->GID(dof) < gmin) gmin = SubMap->GID(dof);
    }
    columns[nC+i] = gmin;
  }
  SubMapC = new Epetra_Map(-1, nC+ndof_set, columns, 0, Comm);
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
  Epetra_Map OwnMapC_temp(-1, ncol_own, columns_own, 0, Comm);
  delete [] columns; delete [] columns_own;
  //
  // determine entries of constraint matrix and substructure dof weights
  //
  get_matrix_diag(A, sub_diag);
  Epetra_Vector Sub_diag(View, *SubMap, sub_diag);
  Epetra_Vector Glob_diag(*OwnMap);
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
  if (MyPID == 0) {
    cout << "weights = " << endl;
    for (i=0; i<nB; i++) cout << weight[i] << endl;
  }
  */
  //
  // solve for substructure coarse interpolation matrices
  //
  ARinvCT = new double[nR*ndof_set];
  double *TEMP = new double[nR*ndof_set];
  if (ndof_set > 0) AR->solve(ndof_set, CT, ARinvCT, TEMP);
  CARinvCT = new double[ndof_set*ndof_set];
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  char TRANSA = 'T'; char TRANSB = 'N'; double ALPHA(1); double BETA(0);
  char UPLO = 'U';
  if (ndof_set > 0)
    EB.GEMM(TRANSA, TRANSB, ndof_set, ndof_set, nR, ALPHA, CT, nR,
	    ARinvCT, nR, BETA, CARinvCT, ndof_set);

  if (MyPID == -1) {
    cout << "CARinvCT = " << endl;
    for (i=0; i<ndof_set; i++) {
      for (j=0; j<ndof_set; j++) cout << CARinvCT[j*ndof_set+i] << " ";
      cout << endl;
    }
  }
  if (ndof_set > 0) {
    EL.POTRF(UPLO, ndof_set, CARinvCT, ndof_set, &INFO);
    assert (INFO == 0);
  }
  delete [] TEMP;
  RHS_cg  = new double[ndof_sub]; 
  SOL_cg  = new double[ndof_sub]; 
  TEMP_cg = new double[ndof_sub];
  lambda_r = new double[ndof_set];
  int ncon = nC + ndof_set;
  double* phi = new double[ndof_sub];
  double* Kc_sub = new double[ncon*ncon];
  PhiB = new double[nB*ncon];
  int* rowbeg = A->get_rowbeg();
  int* colidx = A->get_colidx();
  double* vals = A->get_vals();
  for (i=0; i<ncon; i++) {
    myzero(phi, ndof_sub);
    if (i < nC) {
      phi[dofC[i]] = 1.0;
      myzero(RHS_cg, nR);
      for (j=rowbeg[dofC[i]]; j<rowbeg[dofC[i]+1]; j++) {
	dof = sub_to_R[colidx[j]];
	if (dof >=0) RHS_cg[dof] = vals[j];
      }
      if (nR > 0) AR->solve(1, RHS_cg, SOL_cg, TEMP_cg);
      if (ndof_set > 0) 
	EB.GEMV(TRANSA, nR, ndof_set, ALPHA, CT, nR, SOL_cg, BETA, lambda_r);
    }
    else {
      myzero(lambda_r, ndof_set);
      lambda_r[i-nC] = 1;
      myzero(SOL_cg, nR);
    }
    if (ndof_set > 0) {
      EL.POTRS(UPLO, ndof_set, 1, CARinvCT, ndof_set, lambda_r, ndof_set, 
	       &INFO);
      assert (INFO == 0);
      EB.GEMV(TRANSB, nR, ndof_set, ALPHA, ARinvCT, nR, lambda_r, -1.0, 
	      SOL_cg);
    }
    else EB.SCAL(ndof_set, -1, SOL_cg);
    for (j=0; j<nR; j++) phi[dofR[j]] = SOL_cg[j];
    for (j=0; j<nB; j++) PhiB[nB*i+j] = phi[dofB[j]];
    for (j=0; j<nC; j++) {
      dof = dofC[j]; sum = 0;
      for (k=rowbeg[dof]; k<rowbeg[dof+1]; k++) sum += vals[k]*phi[colidx[k]];
      Kc_sub[ncon*i+j] = sum;
    }
    for (j=0; j<ndof_set; j++) Kc_sub[ncon*i+nC+j] = lambda_r[j];
  }
  delete [] CT; delete [] sub_to_R; delete [] phi;
  if (MyPID == -1) {
    cout << "Kc_sub = " << endl;
    for (i=0; i<ncon; i++) {
      for (j=0; j<ncon; j++) cout << Kc_sub[ncon*j+i] << " ";
      cout << endl;
    }
  }
  //
  // gather coarse stiffness matrix to processor 0
  //
  Epetra_Map Map0(-1, ncol_own, 0, Comm);
  Epetra_IntVector Int_Map0(Map0);
  for (i=0; i<ncol_own; i++) Int_Map0[i] = OwnMapC_temp.GID(i);
  int ngather(0);
  if (MyPID == 0) ngather = OwnMapC_temp.NumGlobalPoints();
  Epetra_Map Map_Gather(-1, ngather, 0, Comm);
  Epetra_IntVector Int_Gather(Map_Gather);
  Epetra_Import Import_Gather(Map_Gather, Map0);
  Int_Gather.Import(Int_Map0, Import_Gather, Insert);
  int *gcdof = new int[ngather];
  for (i=0; i<ngather; i++) gcdof[i] = Int_Gather[i];
  OwnMapC = new Epetra_Map(-1, ngather, gcdof, 0, Comm);
  delete [] gcdof;
  Epetra_CrsMatrix Kc(Copy, *OwnMapC, 0);
  Epetra_CrsMatrix Kc_loc(Copy, *SubMapC, *SubMapC, ncon);
  ImporterC = new Epetra_Import(*SubMapC, *OwnMapC);
  ExporterC = new Epetra_Export(*SubMapC, *OwnMapC);
  double *dvec = new double[ncon];
  int *ivec = new int[ncon];
  for (i=0; i<ncon; i++) ivec[i] = i;
  for (i=0; i<ncon; i++) {
    for (j=0; j<ncon; j++) dvec[j] = Kc_sub[i+j*ncon];
    Kc_loc.InsertMyValues(i, ncon, dvec, ivec);
  }
  delete [] dvec; delete [] ivec; delete [] Kc_sub;
  Kc_loc.FillComplete();
  //  cout << Kc_loc << endl;
  Kc.Export(Kc_loc, *ExporterC, Add);
  Kc.FillComplete(*OwnMapC, *OwnMapC);
  Kc.OptimizeStorage();
  //  cout << Kc << endl;
  //  spmat_datfile(Kc, "Kc.dat", 1);
  //
  // factor coarse stiffness matrix
  //
  if (MyPID == 0) {
    int nrow, nn;
    nrow = Kc.NumGlobalRows();
    cout << "coarse problem dimension = " << nrow << endl;
    if (coarse_solver == 1) {
      rowbeg = new int[nrow+1]; rowbeg[0] =0;
      for (i=0; i<nrow; i++) {
	Kc.ExtractMyRowView(i, nn, vals, colidx);
	rowbeg[i+1] = rowbeg[i] + nn;
      }
      Kc.ExtractMyRowView(0, nn, vals, colidx);
      AKc = new CLAPS_sparse_lu2();
      AKc->factor(vals, rowbeg, colidx, nrow, scale_option);
      delete [] rowbeg;
      SOL_Kc  = new double[nrow];
      TEMP_Kc = new double[nrow];
    }
  }
  //
  // set up maps, importers, and exporters for boundary dofs
  //
  nB_own = 0;
  int* gdofB = new int[nB];
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
  SubMapB = new Epetra_Map(-1, nB, gdofB, 0, Comm);
  OwnMapB = new Epetra_Map(-1, nB_own, gdofB_own, 0, Comm);
  delete [] gdofB; delete [] gdofB_own; delete [] owner_flag;
  ImporterB = new Epetra_Import(*SubMapB, *OwnMapB);
  ExporterB = new Epetra_Export(*SubMapB, *OwnMapB);
}

void CLIP_solver2::allocate_vectors()
{
  int NumVectors(1);
  rSub    = new Epetra_MultiVector(*SubMap,  NumVectors);
  uSub    = new Epetra_MultiVector(*SubMap,  NumVectors);
  rOwn    = new Epetra_MultiVector(*OwnMap,  NumVectors);
  uOwn    = new Epetra_MultiVector(*OwnMap,  NumVectors);
  rSubB   = new Epetra_MultiVector(*SubMapB, NumVectors);
  uSubB   = new Epetra_MultiVector(*SubMapB, NumVectors);
  rOwnB   = new Epetra_MultiVector(*OwnMapB, NumVectors);
  uOwnB   = new Epetra_MultiVector(*OwnMapB, NumVectors);  
  rSubC   = new Epetra_MultiVector(*SubMapC, NumVectors);
  uSubC   = new Epetra_MultiVector(*SubMapC, NumVectors);
  rOwnC   = new Epetra_MultiVector(*OwnMapC, NumVectors);
  uOwnC   = new Epetra_MultiVector(*OwnMapC, NumVectors);
}

double CLIP_solver2::norm2(double a[], int n)
{
  double sum(0), sum_all;
  for (int i=0; i<n; i++) sum += a[i]*a[i];
  Comm.SumAll(&sum, &sum_all, 1);
  return sqrt(sum_all);
}

double CLIP_solver2::dotprod(double a[], double b[], int n)
{
  double sum(0), sum_all;
  for (int i=0; i<n; i++) sum += a[i]*b[i];
  Comm.SumAll(&sum, &sum_all, 1);
  return sum_all;
}

void CLIP_solver2::sum_vectors(double a[], int n, double a_sum[])
{
  Comm.SumAll(a, a_sum, n);
}

int CLIP_solver2::initialize_solve(double u[], double r[])
{
  if (sub_solver <= 1) {
    int i;
    int MyLDA;
    double *rOwn_ptr, *rSub_ptr, *uOwn_ptr, *uSub_ptr;
    rOwn->ExtractView(&rOwn_ptr, &MyLDA);
    rSub->ExtractView(&rSub_ptr, &MyLDA);
    uOwn->ExtractView(&uOwn_ptr, &MyLDA);
    uSub->ExtractView(&uSub_ptr, &MyLDA);
    memcpy(rOwn_ptr, r, rOwn->MyLength()*sizeof(double));
    rSub->Import(*rOwn, *Importer, Insert);
    for (i=0; i<nI; i++) RHS_cg[i] = rSub_ptr[dofI[i]];
    AI->solve(1, RHS_cg, SOL_cg, TEMP_cg);
    myzero(uSub_ptr, ndof_sub);
    for (i=0; i<nI; i++) uSub_ptr[dofI[i]] = SOL_cg[i];
    uOwn->Export(*uSub, *Exporter, Add);
    memcpy(u, uOwn_ptr, uOwn->MyLength()*sizeof(double));
    return 1;
  }
  return 0;
}

void CLIP_solver2::apply_preconditioner(const double r[], double z[])
{
  int i, MyLDA;
  double *rOwn_ptr, *rSub_ptr, *rSubB_ptr, *uOwn_ptr;
  rOwn->ExtractView(&rOwn_ptr, &MyLDA);
  uOwn->ExtractView(&uOwn_ptr, &MyLDA);
  rSub->ExtractView(&rSub_ptr, &MyLDA);
  rSubB->ExtractView(&rSubB_ptr, &MyLDA);
  memcpy(rOwn_ptr, r, rOwn->MyLength()*sizeof(double));
  rSub->Import(*rOwn, *Importer, Insert);
  for (i=0; i<nB; i++) rSubB_ptr[i] = rSub_ptr[dofB[i]];
  uSubB->PutScalar(0);
  //
  // initial static condensation correction
  //
  if (sub_solver >= 1) static_correction1();
  //
  // scale residual
  //
  weight_scale(rSubB);
  //
  // coarse problem correction
  //
  coarse_correction();
  //
  // substructure correction
  //
  substructure_correction();
  //
  // scale preconditioned residual
  //
  weight_scale(uSubB);
  //
  // final static condensation correction
  //
  static_correction2();
  memcpy(z, uOwn_ptr, uOwn->MyLength()*sizeof(double));
}

void CLIP_solver2::static_correction1()
{
  weight_scale(rSubB);
  int i, MyLDA;
  double *rSub_ptr, *rSubB_ptr;
  rSub->ExtractView(&rSub_ptr, &MyLDA);
  rSubB->ExtractView(&rSubB_ptr, &MyLDA);
  for (i=0; i<nI; i++) RHS_cg[i] = rSub_ptr[dofI[i]];
  AI->solve(1, RHS_cg, SOL_cg, TEMP_cg);
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nI; i++) TEMP_cg[dofI[i]] = SOL_cg[i];
  A_times_x_sub(TEMP_cg, RHS_cg, dofB, nB);
  for (i=0; i<nB; i++) rSubB_ptr[i] -= RHS_cg[i];
  rOwnB->PutScalar(0);
  rOwnB->Export(*rSubB, *ExporterB, Add);
  rSubB->Import(*rOwnB, *ImporterB, Insert);
}

void CLIP_solver2::static_correction2()
{
  uOwnB->PutScalar(0);
  uOwnB->Export(*uSubB, *ExporterB, Add);
  uSubB->Import(*uOwnB, *ImporterB, Insert);
  int i, MyLDA;
  double *uSubB_ptr, *uSub_ptr, *rSub_ptr;
  uSubB->ExtractView(&uSubB_ptr, &MyLDA);
  uSub->ExtractView(&uSub_ptr, &MyLDA);
  rSub->ExtractView(&rSub_ptr, &MyLDA);
  myzero(uSub_ptr, ndof_sub);
  for (i=0; i<nB; i++) uSub_ptr[dofB[i]] = uSubB_ptr[i];
  A_times_x_sub(uSub_ptr, RHS_cg, dofBB, nBB);
  for (i=0; i<nBB; i++) rSub_ptr[dofBB[i]] -= RHS_cg[i];
  myzero(RHS_cg, nI);
  for (i=0; i<nI; i++) RHS_cg[i] = rSub_ptr[dofI[i]];
  AI->solve(1, RHS_cg, SOL_cg, TEMP_cg);
  for (i=0; i<nI; i++) uSub_ptr[dofI[i]] = SOL_cg[i];
  uOwn->Export(*uSub, *Exporter, Insert);
}

void CLIP_solver2::coarse_correction()
{
  int MyLDA;
  double *rSubB_ptr, *rSubC_ptr, *rOwnC_ptr, *uSubC_ptr, *uOwnC_ptr;
  double *uSubB_ptr, ALPHA(1), BETA(0);
  char TRANS = 'T';
  Epetra_BLAS EB;
  rSubB->ExtractView(&rSubB_ptr, &MyLDA);
  rSubC->ExtractView(&rSubC_ptr, &MyLDA);
  rOwnC->ExtractView(&rOwnC_ptr, &MyLDA);
  uSubC->ExtractView(&uSubC_ptr, &MyLDA);
  uOwnC->ExtractView(&uOwnC_ptr, &MyLDA);
  uSubB->ExtractView(&uSubB_ptr, &MyLDA);
  EB.GEMV(TRANS, nB, nC+ndof_set, ALPHA, PhiB, nB, rSubB_ptr,
	  BETA, rSubC_ptr);
  rOwnC->PutScalar(0);
  rOwnC->Export(*rSubC, *ExporterC, Add);
  if (MyPID == 0) AKc->solve(1, rOwnC_ptr, uOwnC_ptr, TEMP_Kc);
  uSubC->Import(*uOwnC, *ImporterC, Insert);
  TRANS = 'N';
  EB.GEMV(TRANS, nB, nC+ndof_set, ALPHA, PhiB, nB, uSubC_ptr,
	  BETA, uSubB_ptr);
}

void CLIP_solver2::substructure_correction()
{
  int i, MyLDA, INFO;
  char UPLO('U'), TRANS('T');
  double *rSubB_ptr, *uSubB_ptr, ALPHA(1), BETA(0);
  Epetra_BLAS EB;
  Epetra_LAPACK EL;
  rSubB->ExtractView(&rSubB_ptr, &MyLDA);
  uSubB->ExtractView(&uSubB_ptr, &MyLDA);
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nB; i++) TEMP_cg[dofB[i]] = rSubB_ptr[i];
  myzero(RHS_cg, nR);
  for (i=0; i<nR; i++) RHS_cg[i] = TEMP_cg[dofR[i]];
  AR->solve(1, RHS_cg, SOL_cg, TEMP_cg);
  if (ndof_set > 0) {
    EB.GEMV(TRANS, nR, ndof_set, ALPHA, ARinvCT, nR, RHS_cg,
	    BETA, lambda_r);
    EL.POTRS(UPLO, ndof_set, 1, CARinvCT, ndof_set, lambda_r, ndof_set, &INFO);
    assert (INFO == 0);
    TRANS = 'N'; ALPHA = -1; BETA = 1;
    EB.GEMV(TRANS, nR, ndof_set, ALPHA, ARinvCT, nR, lambda_r,
	    BETA, SOL_cg);
  }
  myzero(TEMP_cg, ndof_sub);
  for (i=0; i<nR; i++) TEMP_cg[dofR[i]] = SOL_cg[i];
  for (i=0; i<nB; i++) uSubB_ptr[i] += TEMP_cg[dofB[i]];
}

void CLIP_solver2::weight_scale(Epetra_MultiVector* Mvec)
{
  int i, MyLDA;
  double* Mvec_ptr;
  Mvec->ExtractView(&Mvec_ptr, &MyLDA);
  assert (Mvec->MyLength() == nB);
  for (i=0; i<nB; i++) Mvec_ptr[i] *= weight[i];
}

void CLIP_solver2::A_times_x_sub(double* x, double* Ax)
{
  int i, j;
  double sum;
  int* rowbeg = A->get_rowbeg();
  int* colidx = A->get_colidx();
  double* vals = A->get_vals();
  int nrow = A->get_nrow();
  for (i=0; i<nrow; i++) {
    sum = 0;
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) {
      sum += vals[j]*x[colidx[j]];
    }
    Ax[i] = sum;
  }
}

void CLIP_solver2::A_times_x_sub(double* x, double* Ax, int* rows, int n)
{
  int i, j, dof;
  double sum;
  int* rowbeg = A->get_rowbeg();
  int* colidx = A->get_colidx();
  double* vals = A->get_vals();
  for (i=0; i<n; i++) {
    dof = rows[i];
    sum = 0;
    for (j=rowbeg[dof]; j<rowbeg[dof+1]; j++) {
      sum += vals[j]*x[colidx[j]];
    }
    Ax[i] = sum;
  }
}

void CLIP_solver2::A_times_x(double* x, double* Ax)
{
    int MyLDA;
    double *rOwn_ptr, *rSub_ptr, *uOwn_ptr, *uSub_ptr;
    rOwn->ExtractView(&rOwn_ptr, &MyLDA);
    rSub->ExtractView(&rSub_ptr, &MyLDA);
    uOwn->ExtractView(&uOwn_ptr, &MyLDA);
    uSub->ExtractView(&uSub_ptr, &MyLDA);
    memcpy(uOwn_ptr, x,  uOwn->MyLength()*sizeof(double));
    uSub->Import(*uOwn, *Importer, Insert);
    A_times_x_sub(uSub_ptr, rSub_ptr);
    rOwn->PutScalar(0);
    rOwn->Export(*rSub, *Exporter, Add);
    memcpy(Ax, rOwn_ptr, rOwn->MyLength()*sizeof(double));
}

void CLIP_solver2::get_matrix_diag(CRS_serial* AA, double* & diag)
{
  int i, j, n;
  n = AA->get_nrow();
  int* rowbeg = AA->get_rowbeg();
  int* colidx = AA->get_colidx();
  double* vals = AA->get_vals();
  diag = new double[n]; myzero(diag, n);
  for (i=0; i<n; i++) {
    for (j=rowbeg[i]; j<rowbeg[i+1]; j++) {
      if (colidx[j] == i) diag[i] = vals[j];
    }
  }
}

void CLIP_solver2::gen_matrix(CRS_serial* AA, int na, int adof[], int* & rowbeg, 
			      int* & colidx, double* &vals)
{
  int i, j, *imap, nnz(0), dof, dof2;
  int n = AA->get_nrow();
  int* rowbegA = AA->get_rowbeg();
  int* colidxA = AA->get_colidx();
  double* valsA = AA->get_vals();
  imap = new int[n]; for (i=0; i<n; i++) imap[i] = -1;
  rowbeg = new int[na+1]; rowbeg[0] = 0;
  for (i=0; i<na; i++) imap[adof[i]] = i;
  for (i=0; i<na; i++) {
    dof = adof[i];
    for (j=rowbegA[dof]; j<rowbegA[dof+1]; j++) {
      dof2 = colidxA[j];
      if (imap[dof2] >= 0) nnz++;
    }
    rowbeg[i+1] = nnz;
  }
  colidx = new int[nnz];
  vals = new double[nnz];
  nnz = 0;
  for (i=0; i<na; i++) {
    dof = adof[i];
    for (j=rowbegA[dof]; j<rowbegA[dof+1]; j++) {
      dof2 = colidxA[j];
      if (imap[dof2] >= 0) {
	colidx[nnz] = imap[dof2];
	vals[nnz] = valsA[j];
	nnz++;
      }
    }
  }
  delete [] imap; 
}

int CLIP_solver2::determine_shell(int node)
{
  int i, flag(0);
  int* loc_dofs = A->get_localdof();
  int* nodebeg = A->get_nodebeg();
  for (i=nodebeg[node]; i<nodebeg[node+1]; i++) {
    if ((loc_dofs[i] >= 4) && (loc_dofs[i] <= 6)) flag = 1;
  }
  return flag;
}

bool CLIP_solver2::find_sub(int dof, int sub)
{
  int i;
  bool answer(false);
  for (i=sub2[dof]; i<sub2[dof+1]; i++) {
    if (sub1[i] == sub) {
      answer = true;
      break;
    }
  }
  return answer;
}

void CLIP_solver2::check_two_nodes(int ii, int* adof, bool* adof_flag,
				   int* anode, bool* anode_flag)
{
  int i, node1, node2, check(0), nadof, nanode;
  int* nodebeg = A->get_nodebeg();
  double* x = A->get_xcoord();
  double* y = A->get_ycoord();
  int* aa = &dset1[dset2[ii]];
  int  nn = dset2[ii+1] - dset2[ii];
  get_adof(ii, aa, nn, adof, nadof, adof_flag);
  get_anode(adof, nadof, anode, anode_flag, nanode);
  if (nanode <= 2) {
    for (i=0; i<nanode; i++) set_corner_flag(anode[i]);
    return;
  }
  for (i=0; i<nanode; i++)
    if (corner_flag[nodebeg[anode[i]]] == 1) check++;
  //  if (MyPID == 0) cout << "check = " << check << endl;
  if (check >= 2) return;
  if (check == 1) {
    node1 = find_node1(anode, nanode);
    find_farthest1(anode, nanode, node1, node2);
    set_corner_flag(node2);
    return;
  }
  if (check == 0) {
    find_farthest2(anode, nanode, node1, node2);
    set_corner_flag(node1);
    set_corner_flag(node2);
  }
}
    
void CLIP_solver2::set_corner_flag(int node)
{
  int* nodebeg = A->get_nodebeg();
  for (int i=nodebeg[node]; i<nodebeg[node+1]; i++) corner_flag[i] = 1;
}

void CLIP_solver2::check_three_nodes(int ii, int* adof, bool* adof_flag,
				    int* anode, bool* anode_flag)
{
  int i, node1, node2, node3, check(0), nadof, nanode;
  int* nodebeg = A->get_nodebeg();
  int* aa = &dset1[dset2[ii]];
  int  nn = dset2[ii+1] - dset2[ii];
  get_adof(ii, aa, nn, adof, nadof, adof_flag);
  get_anode(adof, nadof, anode, anode_flag, nanode);
  //  if (MyPID == 0) cout << "nanode = " << nanode << endl;
  if (nanode <= 3) {
    for (i=0; i<nanode; i++) set_corner_flag(anode[i]);
    return;
  }
  for (i=0; i<nanode; i++)
    if (corner_flag[nodebeg[anode[i]]] == 1) check++;
  if (check >= 3) return;
  if (check == 2) {
    find_node1_node2(anode, nanode, node1, node2);
    node3 = find_max_area(anode, nanode, node1, node2);
    set_corner_flag(node3);
    return;
  }
  if (check == 1) {
    node1 = find_node1(anode, nanode);
    find_farthest1(anode, nanode, node1, node2);
    node3 = find_max_area(anode, nanode, node1, node2);
    set_corner_flag(node2); set_corner_flag(node3);
    return;
  }
  if (check == 0) {
    find_farthest2(anode, nanode, node1, node2);
    node3 = find_max_area(anode, nanode, node1, node2);
    set_corner_flag(node1); set_corner_flag(node2); set_corner_flag(node3);
    return;
  }
}

void CLIP_solver2::get_adof(int ii, int* aa, int nn, int* adof, int & nadof, 
			    bool* adof_flag, int tflag)
{
  int i, j, dof, dof2, flag;
  int *rowbeg = A->get_rowbeg();
  int *colidx = A->get_colidx();
  nadof = 0;
  for (i=0; i<nn; i++) {
    dof = aa[i];
    adof[nadof] = dof;
    adof_flag[dof] = true;
    nadof++;
  }
  int nadof_orig = nadof;
  for (i=0; i<nadof_orig; i++) {
    dof = adof[i];
    for (j=rowbeg[dof]; j<rowbeg[dof+1]; j++) {
      dof2 = colidx[j];
      if (adof_flag[dof2] == false) {
	flag = 0;
	if (tflag == 1) flag = 1;
	else if (contains(ii, dof2) == true) flag = 1;
	if (flag == 1) {
	  adof[nadof] = dof2;
	  adof_flag[dof2] = true;
	  nadof++;
	}
      }
    }
  }
  for (i=0; i<nadof; i++) adof_flag[adof[i]] = false;
}

bool CLIP_solver2::contains(int ii, int dof2)
{
  int i, dof, count(0);
  for (i=sub2[dof2]; i<sub2[dof2+1]; i++) sub_flag[sub1[i]] = true;
  dof = dset1[dset2[ii]];
  for (i=sub2[dof]; i<sub2[dof+1]; i++) 
    if (sub_flag[sub1[i]] == true) count++;
  for (i=sub2[dof2]; i<sub2[dof2+1]; i++) sub_flag[sub1[i]] = false;
  if (count == (sub2[dof+1]-sub2[dof])) return true;
  else return false;
}

void CLIP_solver2::get_anode(int* adof, int nadof, int* anode, 
			     bool* anode_flag, int & nanode)
{
  int i, node;
  nanode = 0;
  for (i=0; i<nadof; i++) {
    node = dof2node[adof[i]];
    if (anode_flag[node] == false) {
      anode[nanode] = node;
      anode_flag[node] = true;
      nanode++;
    }
  }
  for (i=0; i<nanode; i++) anode_flag[anode[i]] = false;
}

int CLIP_solver2::find_node1(int* anode, int nanode)
{
  int i, check(0), node, node1;
  int* nodebeg = A->get_nodebeg();
  for (i=0; i<nanode; i++) {
    node = anode[i];
    if (corner_flag[nodebeg[node]] == 1) {
      node1 = node;
      check++;
    }
  }
  assert(check == 1);
  return node1;
}

void CLIP_solver2::find_farthest1(int* anode, int nanode, int node1, int & node2)
{
  int i, node;
  int* nodebeg = A->get_nodebeg();
  double* x = A->get_xcoord();
  double* y = A->get_ycoord();
  double* z = A->get_zcoord();
  double dist2, max_dist2(-1), dx, dy, dz;
  node2 = -1;
  for (i=0; i<nanode; i++) {
    node = anode[i];
    if ((corner_flag[nodebeg[node]] == 0) && (node != node1)) {
      dx = x[node] - x[node1];
      dy = y[node] - y[node1];
      dz = z[node] - z[node1];
      dist2 = dx*dx + dy*dy + dz*dz;
      if (dist2 > max_dist2) {
	node2 = node;
	max_dist2 = dist2;
      }
    }
  }
  assert(node2 != -1);
}

void CLIP_solver2::find_farthest2(int* anode, int nanode, int & node1, int & node2)
{
  int i, j, nodec1, nodec2;
  int* nodebeg = A->get_nodebeg();
  double* x = A->get_xcoord();
  double* y = A->get_ycoord();
  double* z = A->get_zcoord();
  double dist2, max_dist2(-1), dx, dy, dz;
  node1 = -1;
  for (i=0; i<nanode; i++) {
    nodec1 = anode[i];
    if (corner_flag[nodebeg[nodec1]] == 0) {
      for (j=1; j<nanode; j++) {
	nodec2 = anode[j];
	if (corner_flag[nodebeg[nodec2]] == 0) {
	  dx = x[nodec2] - x[nodec1];
	  dy = y[nodec2] - y[nodec1];
	  dz = z[nodec2] - z[nodec1];
	  dist2 = dx*dx + dy*dy + dz*dz;
	  if (dist2 > max_dist2) {
	    node1 = nodec1;
	    node2 = nodec2;
	    max_dist2 = dist2;
	  }
	}
      }
    }
  }
  assert(node1 != -1);
}

void CLIP_solver2::find_node1_node2(int* anode, int nanode, int & node1, int & node2)
{
  int i, check(0), node;
  int* nodebeg = A->get_nodebeg();
  for (i=0; i<nanode; i++) {
    node = anode[i];
    if (corner_flag[nodebeg[node]] == 1) {
      check++;
      if (check == 1) node1 = node;
      if (check == 2) node2 = node;
    }
  }
  assert(check == 2);
}

int CLIP_solver2::find_max_area(int* anode, int nanode, int node1, int node2)
{
  int i, node, node3(-1);
  int* nodebeg = A->get_nodebeg();;
  double* x = A->get_xcoord();
  double* y = A->get_ycoord();
  double* z = A->get_zcoord();
  double area2, max_area2(-1), dx1, dy1, dz1, dx2, dy2, dz2, cx, cy, cz;
  dx1 = x[node2] - x[node1];
  dy1 = y[node2] - y[node1];
  dz1 = z[node2] - z[node1];
  for (i=0; i<nanode; i++) {
    node = anode[i];
    if (corner_flag[nodebeg[node]] == 0) {
      dx2 = x[node] - x[node1];
      dy2 = y[node] - y[node1];
      dz2 = z[node] - z[node1];
      cx = dy1*dz2 - dz1*dy2;
      cy = dz1*dx2 - dx1*dz2;
      cz = dx1*dy2 - dy1*dx2;
      area2 = cx*cx + cy*cy + cz*cz;
      if (area2 > max_area2) {
	max_area2 = area2;
	node3 = node;
      }
    }
  }
  assert(node3 != -1);
  return node3;
}

void CLIP_solver2::determine_components(int A1[], int A2[], int N, 
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

void CLIP_solver2::zero_pointers()
{
  AKc = 0; SOL_Kc = 0; TEMP_Kc = 0;
  amg_matR = 0; amg_A1R = 0; amg_A2R = 0; amg_xR = 0;
  amg_yR = 0; amg_zR = 0; amg_nodebegR = 0; amg_local_dofR = 0;
  amg_matI = 0; amg_A1I = 0; amg_A2I = 0; amg_xI = 0;
  amg_yI = 0; amg_zI = 0; amg_nodebegI = 0; amg_local_dofI = 0;
  amg_matC = 0; amg_A1C = 0; amg_A2C = 0; amg_xC = 0;
  amg_yC = 0; amg_zC = 0; amg_nodebegC = 0; amg_local_dofC = 0;
}


void CLIP_solver2::spmat_datfile(const Epetra_CrsMatrix & AA, char fname[], 
				int opt)
{
  int iproc, i, j, NumEntries, *Indices, grow, gcol;
  double *Values;
  ofstream ffout;
  for (iproc=0; iproc<AA.Comm().NumProc(); iproc++) {
    if (AA.Comm().MyPID() == iproc) {
      if (AA.Comm().MyPID() == 0) ffout.open(fname);
      else ffout.open(fname, ofstream::out | ofstream::app);
      for (i=0; i<AA.NumMyRows(); i++) { 
	AA.ExtractMyRowView(i, NumEntries, Values, Indices);
	for (j=0; j<NumEntries; j++) {
	  if (opt == 1)
	    ffout << i+1 << " " << Indices[j]+1 << setw(22) << setprecision(15)
		  << Values[j] << endl;
	  if (opt == 2) {
	    grow = AA.GRID(i); gcol = AA.GCID(Indices[j]);
	    ffout << grow+1 << " " << gcol+1 << setw(22) 
		  << setprecision(15) << Values[j] << endl;
	  }
	}
      }
      ffout.close();
    }
    AA.Comm().Barrier();
  }
}

void CLIP_solver2::gen_amg_data(int n, int* dofa, double* & amg_mat, int* & amg_A1,
				int* & amg_A2, double* & amg_x, double* & amg_y,
				double* & amg_z, int* & amg_nodebeg,
				int* & amg_local_dof, int & amg_nnode)
{
  int i, j, k, node, prev_node(-1);
  int* local_dof = A->get_localdof();
  double *x = A->get_xcoord();
  double *y = A->get_ycoord();
  double *z = A->get_zcoord();
  int nnode = A->get_nnode();
  int* rowbeg = A->get_rowbeg();
  int* colidx = A->get_colidx();
  double* vals = A->get_vals();
  amg_local_dof = new int[n];
  amg_nnode = 0;
  for (i=0; i<n; i++) {
    amg_local_dof[i] = local_dof[dofa[i]];
    node = dof2node[dofa[i]];
    if (node != prev_node) {
      amg_nnode++;
      prev_node = node;
    }
  }
  amg_nodebeg = new int[amg_nnode+1];
  amg_x = new double[amg_nnode];
  amg_y = new double[amg_nnode];
  amg_z = new double[amg_nnode];
  amg_nnode = 0; 
  prev_node = -1;
  for (i=0; i<n; i++) {
    node = dof2node[dofa[i]];
    if (node != prev_node) {
      amg_nodebeg[amg_nnode] = i;
      amg_x[amg_nnode] = x[node];
      amg_y[amg_nnode] = y[node];
      amg_z[amg_nnode] = z[node];
      amg_nnode++;
      prev_node = node;
    }
  }
  amg_nodebeg[amg_nnode] = n;
  int* imap = new int[nnode];
  for (i=0; i<nnode; i++) imap[i] = -1;
  amg_nnode = 0;
  prev_node = -1;
  for (i=0; i<n; i++) {
    node = dof2node[dofa[i]];
    if (node != prev_node) {
      imap[node] = amg_nnode;
      amg_nnode++;
      prev_node = node;
    }
  }
  amg_A2 = new int[amg_nnode+1]; amg_A2[0] = 0;
  int* anode = new int[amg_nnode];
  bool* flag = new bool[amg_nnode];
  for (i=0; i<amg_nnode; i++) flag[i] = false;
  int nanode, dof1, dof2, amg_node;
  //
  // pass 1 to determine amg_A2
  //
  for (i=0; i<amg_nnode; i++) {
    nanode = 0;
    for (j=amg_nodebeg[i]; j<amg_nodebeg[i+1]; j++) {
      dof1 = dofa[j];
      for (k=rowbeg[dof1]; k<rowbeg[dof1+1]; k++) {
	dof2 = colidx[k];
	node = dof2node[dof2];
	amg_node = imap[node];
	if (amg_node != -1) {
	  if (flag[amg_node] == false) {
	    anode[nanode] = amg_node;
	    flag[amg_node] = true;
	    nanode++;
	  }
	}
      }
    }
    amg_A2[i+1] = amg_A2[i] + nanode;
    for (j=0; j<nanode; j++) flag[anode[j]] = false;
  }
  //
  // pass 2 to determine amg_A1 and amg_A1_beg
  //
  amg_A1 = new int[amg_A2[amg_nnode]];
  int *amg_A1_beg = new int[amg_A2[amg_nnode]+1]; amg_A1_beg[0] = 0;
  int nnz(0), ndofpn1, ndofpn2;
  for (i=0; i<amg_nnode; i++) {
    ndofpn1 = amg_nodebeg[i+1] - amg_nodebeg[i];
    nanode = 0;
    for (j=amg_nodebeg[i]; j<amg_nodebeg[i+1]; j++) {
      dof1 = dofa[j];
      for (k=rowbeg[dof1]; k<rowbeg[dof1+1]; k++) {
	dof2 = colidx[k];
	node = dof2node[dof2];
	amg_node = imap[node];
	if (amg_node != -1) {
	  if (flag[amg_node] == false) {
	    anode[nanode] = amg_node;
	    flag[amg_node] = true;
	    amg_A1[amg_A2[i]+nanode] = amg_node;
	    ndofpn2 = amg_nodebeg[amg_node+1] - amg_nodebeg[amg_node];
	    nnz += ndofpn1*ndofpn2;
	    amg_A1_beg[amg_A2[i]+nanode+1] = nnz;
	    nanode++;
	  }
	}
      }
    }
    for (j=0; j<nanode; j++) flag[anode[j]] = false;
  }
  //
  // pass 3 to determine amg_mat
  //
  if (MyPID == -1) {
    cout << "amg_A1 = " << endl;
    for (i=0; i<amg_nnode; i++) {
      for (j=amg_A2[i]; j<amg_A2[i+1]; j++) cout << amg_A1[j] << " ";
      cout << endl;
    }
    //    cout << "amg_A1_beg = " << endl;
    //    for (i=0; i<amg_nnode; i++) {
    //      for (j=amg_A2[i]; j<amg_A2[i+1]; j++) cout << amg_A1_beg[j] << " ";
    //      cout << endl;
    //    }
    //    cout << amg_A1_beg[amg_A2[amg_nnode]] << endl;
  }
  int ldof, m, loc_col, jj;
  amg_mat = new double[nnz]; myzero(amg_mat, nnz);
  for (i=0; i<amg_nnode; i++) anode[i] = -1;
  for (i=0; i<amg_nnode; i++) {
    for (j=amg_A2[i]; j<amg_A2[i+1]; j++) anode[amg_A1[j]] = j - amg_A2[i];
    for (j=amg_nodebeg[i]; j<amg_nodebeg[i+1]; j++) {
      jj = j - amg_nodebeg[i];
      dof1 = dofa[j];
      for (k=rowbeg[dof1]; k<rowbeg[dof1+1]; k++) {
	dof2 = colidx[k];
	node = dof2node[dof2];
	amg_node = imap[node];
	if (amg_node != -1) {
	  loc_col = anode[amg_node];
	  assert (loc_col != -1);
	  ldof = -1;
	  for (m=amg_nodebeg[amg_node]; m<amg_nodebeg[amg_node+1]; m++) {
	    if (amg_local_dof[m] == local_dof[dof2]) {
	      ldof = m - amg_nodebeg[amg_node];
	      break;
	    }
	  }
	  if (ldof != -1) {
	    amg_mat[amg_A1_beg[amg_A2[i] + loc_col] + ndofpn1*ldof + jj] = vals[k];
	  }
	}
      }
    }
    for (j=amg_A2[i]; j<amg_A2[i+1]; j++) anode[amg_A1[j]] = -1;
  }
  if (MyPID == -1) {
    cout << "amg_mat = " << endl;
    for (i=0; i<amg_nnode; i++) {
      ndofpn1 = amg_nodebeg[i+1] - amg_nodebeg[i];
      for (m=0; m<ndofpn1; m++) {
	for (j=amg_A2[i]; j<amg_A2[i+1]; j++) {
	  node = amg_A1[j];
	  ndofpn2 = amg_nodebeg[node+1] - amg_nodebeg[node];
	  for (int mm=0; mm<ndofpn2; mm++) {
	    cout << amg_mat[amg_A1_beg[j]+m+ndofpn1*mm] << " ";
	  }
	}
	cout << endl;
      }
    }
  }
  delete [] imap; delete [] anode; delete [] flag; delete [] amg_A1_beg;
}
