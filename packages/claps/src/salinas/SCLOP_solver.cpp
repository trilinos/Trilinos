#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "SCLOP_solver.hpp"
#include "myzero.hpp"

SCLOP_solver::SCLOP_solver(
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
  AStandard = 0; ConStandard = 0; CStandard = 0; SubMap = 0; StandardMap = 0;
  LDStandard = 0; GNStandard = 0; H1 = 0; uStand = 0; uLocal = 0;
  fStand = 0; uSub = 0; uvec = 0; fvec = 0; subvec = 0; locvec = 0;
  ImporterSt2Sub = 0; ExporterSub2St = 0; MpcLocalMap = 0;
  ImporterStLam2Loc = 0; CS = 0;
}

SCLOP_solver::~SCLOP_solver()
{
  delete Comm; delete AStandard; delete ConStandard; delete CStandard;
  delete LDStandard; delete GNStandard; delete uStand; delete uLocal;
  delete fStand; delete uSub; delete [] uvec; delete [] fvec;
  delete [] subvec; delete [] locvec; delete SubMap; delete StandardMap; 
  delete ImporterSt2Sub; delete ExporterSub2St; delete MpcLocalMap;
  delete ImporterStLam2Loc; delete CS;
}

void SCLOP_solver::construct_H1()
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
      if (n > 1) sort(&H1[H2[i]], &H1[H2[i+1]]);
    }
    //
    // sort global node numbers
    //
    int *gnn_sort = new int[nnode];
    for (i=0; i<nnode; i++) gnn_sort[i] = gnn[i];
    sort(&gnn_sort[0], &gnn_sort[nnode]);
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

void SCLOP_solver::determine_ownership()
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
  int *count1 = new int[ndof]; myzero(count1, ndof);
  nmpc_loc = 0;
  for (i=0; i<NumberMpc; i++) {
    if (MpcVector[i].NumEntries() > 0) nmpc_loc++;
    for (j=0; j<MpcVector[i].NumEntries(); j++) {
      node = MpcVector[i].LocalId(j);
      dof  = MpcVector[i].NodalDof(j);
      for (k=nodebeg[node]; k<nodebeg[node+1]; k++) {
	if (nodedof[k] == dof) {
	  count1[k]++;
	  break;
	}
      }
    }
  }
  int *con2 = new int[ndof+1]; con2[0] = 0;
  for (i=0; i<ndof; i++) con2[i+1] = con2[i] + count1[i];
  myzero(count1, ndof);
  int *con1 = new int[con2[ndof]];
  double *coef = new double[con2[ndof]];
  int *mpc_loc = new int[nmpc_loc]; 
  locvec = new double[nmpc_loc]; nmpc_loc = 0;
  for (i=0; i<NumberMpc; i++) {
    if (MpcVector[i].NumEntries() > 0) {
      mpc_loc[nmpc_loc] = i;
      nmpc_loc++;
    }
    for (j=0; j<MpcVector[i].NumEntries(); j++) {
      node = MpcVector[i].LocalId(j);
      dof  = MpcVector[i].NodalDof(j);
      for (k=nodebeg[node]; k<nodebeg[node+1]; k++) {
	if (nodedof[k] == dof) {
	  con1[con2[k] + count1[k]] = i;
	  coef[con2[k] + count1[k]] = MpcVector[i].Coef(j);
	  count1[k]++;
	  break;
	}
      }
    }
  }
  Epetra_Map RowMapCon(NumberMpc, 0, *Comm);
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
  double *xyzdof = new double[3*ndof_mine];
  int *localdofs = new int[ndof_mine];
  int *globalnodes = new int[ndof_mine];
  int ndof_mine_keep = ndof_mine;
  ndof_mine = 0;
  int jlow = 0;
  for (i=0; i<nnode; i++) {
    m = maxdofnode*i;
    for (j=jlow; j<nodebeg[i+1]; j++) {
      dof = nodedof[j];
      if (iown[m+dof] == true) {
	xyzdof[ndof_mine                   ] = x[i];
	xyzdof[ndof_mine +   ndof_mine_keep] = y[i];
	xyzdof[ndof_mine + 2*ndof_mine_keep] = z[i];
	localdofs[ ndof_mine] = dof + 1;
	globalnodes[ndof_mine] = gnn[i];
	nodedof[ndof_mine] = dof;
	ndof_mine++;
      }
    }
    jlow = nodebeg[i+1];
    nodebeg[i+1] = ndof_mine;
  }
  delete [] iown; delete [] nodebeg; delete [] nodedof; 

  StandardMap = new Epetra_Map(-1, ndof_mine, sub_gdof, 0, *Comm);
  ConStandard = new Epetra_CrsMatrix(Copy, *StandardMap, 0);
  ConSubdomain.FillComplete(RowMapCon, *StandardMap);
  Epetra_Export Exporter(*SubMap, *StandardMap);
  ConStandard->Export(ConSubdomain, Exporter, Add);
  //  ConStandard->Export(ConSubdomain, Exporter, Insert);
  ConStandard->FillComplete(RowMapCon, *StandardMap);
  //  cout << *ConStandard << endl;

  delete [] count1; delete [] con1; delete [] con2; delete [] coef;
  CStandard = new Epetra_MultiVector(Copy, *StandardMap, xyzdof, ndof_mine, 3);
  LDStandard = new Epetra_IntVector(Copy, *StandardMap, localdofs);
  GNStandard = new Epetra_IntVector(Copy, *StandardMap, globalnodes);
  uStand = new Epetra_Vector(View, *StandardMap, uvec);
  fStand = new Epetra_Vector(View, *StandardMap, fvec);
  uLocal = new Epetra_Vector(View, *MpcLocalMap, locvec);
  delete [] sub_gdof; delete [] xyzdof; delete [] localdofs; 
  delete [] globalnodes;
}

void SCLOP_solver::construct_K_base()
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
  Epetra_CrsMatrix SubdomainMatrix(View, *SubMap, *SubMap, count);
  for (i=0; i<neq; i++) {
    ibeg = rowbeg_base[i];
    ierr = SubdomainMatrix.InsertMyValues(i, count[i], &K_base[ibeg], 
					  &colidx_base[ibeg]);
  }
  SubdomainMatrix.FillComplete();
  AStandard = new Epetra_CrsMatrix(Copy, *StandardMap, 0);
  Epetra_Export Exporter(*SubMap, *StandardMap);
  AStandard->Export(SubdomainMatrix, Exporter, Add);
  AStandard->FillComplete();

  int NumEntries, *Indices;
  double min_diag(1e20), max_diag(0), min_diag_all, max_diag_all, *Values;
  for (i=0; i<AStandard->NumMyRows(); i++) {
    AStandard->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      if (Indices[j] == i) {
	if (Values[j] < min_diag) min_diag = Values[j];
	if (Values[j] > max_diag) max_diag = Values[j];
      }
    }
  }
  Comm->MinAll(&min_diag, &min_diag_all, 1);
  Comm->MaxAll(&max_diag, &max_diag_all, 1);
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
    sprintf(filename,"%s%d","CLOP_base", MyPID);
    CRD_utils::spmat_datfile(neq ,rowbeg_base, colidx_base, K_base, filename);
    EPmat_datfile(AStandard, filename);
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

void SCLOP_solver::CLOP_solver_init(int overlap, double solver_tol,
   int maxiter, int atype, int ndim, int local_solver, int max_orthog,
   int prt_debug, int prt_summary)
{
  double clop_params[20];
  clop_params[0] = double(overlap);
  clop_params[1] = double(solver_tol);
  clop_params[2] = double(maxiter);
  clop_params[3] = double(max_orthog);
  clop_params[4] = double(atype);
  clop_params[5] = double(ndim);
  clop_params[6] = double(local_solver);
  clop_params[7] = double(prt_debug);
  clop_params[8] = double(prt_summary);
  CS = new CLOP_solver(AStandard, LDStandard, CStandard, SubMap, ConStandard, 
 	               GNStandard, clop_params);
}

void SCLOP_solver::solve(double f[], double u[], int & number_iterations, 
			 int & SCLOP_status)
{
  int i;
  for (i=0; i<ndof; i++) subvec[i] = f[i];
  fStand->Export(*uSub, *ExporterSub2St, Add);
  CS->solve(uStand, fStand, number_iterations, SCLOP_status);
  uSub->Import(*uStand, *ImporterSt2Sub, Insert);
  for (i=0; i<ndof; i++) u[i] = subvec[i];
}

void SCLOP_solver::MpcForces( double *cvals)
{
  int i;
  CS->mpcforces(uLocal, ImporterStLam2Loc);
  for (i=0; i<nmpc_loc; i++) cvals[i] = locvec[i];
}

void SCLOP_solver::EPmat_datfile(Epetra_CrsMatrix* A, char fname[])
{
  int i, j, NumEntries, *Indices;
  double *Values;
  ofstream fout;
  sprintf(fname, "%s.epetra", fname);
  fout.open(fname);
  for (i=0; i<A->NumMyRows(); i++) {
    A->ExtractMyRowView(i, NumEntries, Values, Indices);
    for (j=0; j<NumEntries; j++) {
      fout << A->GRID(i)+1 << " " << A->GCID(Indices[j])+1 << setw(22) 
	   << setprecision(15) << Values[j] << endl;
    }
  }
  fout.close();
}
