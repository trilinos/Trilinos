#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <math.h>
#include <assert.h>
#include "CLOP_graph.hpp"
#include "myzero.hpp"
#include "CRD_utils.hpp"

CLOP_graph::CLOP_graph(     
     const Epetra_CrsMatrix* A_,   // stiffness matrix
     const Epetra_IntVector* ND_,  // nodes for dofs
     const int overlap_,           // overlap
     const int partition_option_,  // partitioning option
     const int atype_,             // analysis type
     const int ndim_,              // spatial dimension
     int iwork1[],                 // integer work array
     int iwork2[],                 // integer work array
     int* & dofpart1,              // dofpart1[dofpart2[i]:dofpart2[i+1]-1] =
     int* & dofpart2,              //   dofs in overlapping subdomain i
     int & npart)                  // number of subdomains for processor
  : A(A_), ND(ND_), overlap(overlap_), partition_option(partition_option_),
    atype(atype_), ndim(ndim_), count1(iwork1), imap(iwork2)
{
  //
  // determine target number of dofs for each subdomain
  //
  int min_ndof_sub(10), ncdof_sub, nnode, ncomp, ndof_max, NumProc, MyPID;
  double scale_factor(1.5);
  if (atype == 1) ncdof_sub = 1;
  if (atype == 2) {
    if (ndim == 2) ncdof_sub = 3;
    if (ndim == 3) ncdof_sub = 6;
  }
  if (atype == 3) ncdof_sub = 3;
  ndof = A->NumMyRows();
  A->Comm().MaxAll(&ndof, &ndof_max, 1);
  MyPID = A->Comm().MyPID();
  NumProc = A->Comm().NumProc();
  ndof_target = int(sqrt(1.0*ndof_max*NumProc*ncdof_sub));
  ndof_target = ndof;
  if (ndof_target < min_ndof_sub) ndof_target = min_ndof_sub;
  //  cout << "MyPID, ndof, ndof_target = " << MyPID << " " << ndof << " "
  //       << ndof_target << endl;
  //
  // construct node graph and determine dofs for each node
  //   nnode = number of nodes
  //   node_con1[node_con2[i]:node_con2[i+1]-1] = nodes connected to node i
  //   dof_for_node1[dof_for_node2[i]:dof_for_node2[i+1]-1] = dofs for node i
  //
  int *node_con1, *node_con2, *dof_for_node1, *dof_for_node2;
  node_con1 = node_con2 = dof_for_node1 = dof_for_node2 = 0;
  construct_node_graph(node_con1, node_con2, dof_for_node1, dof_for_node2,
		       nnode);
  //
  // determine components
  //   ncomp = number of components for original node graph
  //   comp1[comp2[i]:comp2[i+1]-1] = nodes in component i
  //
  int *comp1, *comp2;
  comp1 = comp2 = 0;
  determine_components(node_con1, node_con2, nnode, comp1, comp2, ncomp);
  //  cout << "MyPID, ncomp = " << MyPID << " " << ncomp << endl;
  //
  // determine subdomains
  //   nsub = number of subdomains
  //   nsub1[nsub2[i]:nsub2[i+1]-1] = nodes in subdomain i
  //
  int *nsub1, *nsub2;
  nsub1 = nsub2 = 0;
  determine_subdomains(comp1, comp2, ncomp, node_con1, node_con2, 
		       dof_for_node2, nnode, nsub1, nsub2);
  delete [] comp1; delete [] comp2;  
  npart = nsub;
  //  cout << "MyPID, nsub = " << MyPID << " " << nsub << endl;
  //
  // determine overlapping subdomains
  //   nosub1[nosub2[i]:nosub2[i+1]-1] = nodes in overlapping subdomain i
  //
  int *nosub1, *nosub2;
  nosub1 = nosub2 = 0;
  determine_overlap(nsub1, nsub2, node_con1, node_con2, nnode,
		    nosub1, nosub2);
  delete [] node_con1; delete [] node_con2;  
  delete [] nsub1; delete [] nsub2; 
  //
  // determine dofs in overlapping subdomains
  //
  determine_dof_overlap(nosub1, nosub2, nnode, dof_for_node1, dof_for_node2,
			dofpart1, dofpart2);
  delete [] nosub1; delete [] nosub2; delete [] dof_for_node1;
  delete [] dof_for_node2;
}

CLOP_graph::~CLOP_graph()
{
}

void CLOP_graph::construct_node_graph(int* & node_con1, int* & node_con2,
           int* & dof_for_node1, int* & dof_for_node2, int & nnode)
{
  int i, j, k, NumEntries, *Indices, node, nanode, dof;
  double *Values;
  assert (ND->MyLength() == ndof);
  int *node_vec = new int[ndof]; ND->ExtractCopy(node_vec);
  int *dof_map  = new int[ndof]; ND->ExtractCopy(dof_map);
  CRD_utils::sort_and_cull(node_vec, ndof, nnode);
  myzero(count1, nnode);
  for (i=0; i<ndof; i++) {
    node = CRD_utils::find_index(node_vec, nnode, dof_map[i]);
    assert (node != -1);
    dof_map[i] = node;
    count1[node]++;
  }
  delete [] node_vec;
  dof_for_node2 = new int[nnode+1]; dof_for_node2[0] = 0;
  for (i=0; i<nnode; i++) {
    dof_for_node2[i+1] = dof_for_node2[i] + count1[i];
  }
  dof_for_node1 = new int[dof_for_node2[nnode]];
  myzero(count1, nnode);
  for (i=0; i<ndof; i++) {
    node = dof_map[i];
    dof_for_node1[dof_for_node2[node] + count1[node]] = i;
    count1[node]++;
  }
  myzero(count1, nnode);
  node_con2 = new int[nnode+1]; node_con2[0] = 0;
  for (i=0; i<nnode; i++) {
    nanode = 0;
    for (j=dof_for_node2[i]; j<dof_for_node2[i+1]; j++) {
      dof = dof_for_node1[j];
      A->ExtractMyRowView(dof, NumEntries, Values, Indices);
      for (k=0; k<NumEntries; k++) {
	node = dof_map[Indices[k]];
	if (count1[node] == 0) {
	  count1[node] = 1;
	  imap[nanode] = node;
	  nanode++;
	}
      }
    }
    node_con2[i+1] = node_con2[i] + nanode;
    for (j=0; j<nanode; j++) count1[imap[j]] = 0;
  }
  node_con1 = new int[node_con2[nnode]];
  for (i=0; i<nnode; i++) {
    nanode = 0;
    for (j=dof_for_node2[i]; j<dof_for_node2[i+1]; j++) {
      dof = dof_for_node1[j];
      A->ExtractMyRowView(dof, NumEntries, Values, Indices);
      for (k=0; k<NumEntries; k++) {
	node = dof_map[Indices[k]];
	if (count1[node] == 0) {
	  count1[node] = 1;
	  imap[nanode] = node;
	  nanode++;
	}
      }
    }
    for (j=0; j<nanode; j++) {
      count1[imap[j]] = 0;
      node_con1[node_con2[i] + j] = imap[j];
    }
  }
  delete [] dof_map;
}

void CLOP_graph::determine_components(int A1[], int A2[], int N, 
                      int* & comp1, int* & comp2, int & ncomp)
{
  int i, ic;
  CRD_utils::Graph_class Graph(N, A1, A2);
  int *component = new int[N];
  Graph.Components(component, ncomp);
  comp1 = new int[N];
  comp2 = new int[ncomp+1]; comp2[0] = 0;
  int *count_comp = new int[ncomp]; myzero(count_comp, ncomp);
  for (i=0; i<N; i++) count_comp[component[i]]++;
  for (i=0; i<ncomp; i++) comp2[i+1] = comp2[i] + count_comp[i];
  myzero(count_comp, ncomp);
  for (i=0; i<N; i++) {
    ic = component[i];
    comp1[comp2[ic] + count_comp[ic]] = i;
    count_comp[ic]++;
  }
  delete [] component; delete [] count_comp;
}

void CLOP_graph::determine_subdomains(int comp1[], int comp2[], int ncomp, 
	     int node_con1[], int node_con2[], int dof_for_node2[], 
	     int nnode, int* & nsub1, int* & nsub2)
{
  int i, j, k, ii, node, nnode_comp, max_nnode_comp(0), nnz, max_nnz(0);
  int node2, *vwgt, *adjwgt, wgtflag(0), numflag(0), nparts, options[5];
  int edgecut, max_nparts(0);
  vwgt = adjwgt = 0; options[0] = 0;
  int *ndof_comp = new int[ncomp];
  memset(imap, -1, nnode*sizeof(int));
  nsub = 0;
  for (i=0; i<ncomp; i++) {
    nnode_comp = comp2[i+1] - comp2[i];
    if (nnode_comp > max_nnode_comp) max_nnode_comp = nnode_comp;
    ndof_comp[i] = 0;
    for (j=comp2[i]; j<comp2[i+1]; j++) {
      node = comp1[j];
      ndof_comp[i] += dof_for_node2[node+1] - dof_for_node2[node];
      imap[node] = j-comp2[i];
    }
    nparts = ndof_comp[i]/ndof_target;
    if (nparts > max_nparts) max_nparts = nparts;
    if (nparts > 1) {
      nsub += nparts;
      nnz = 0;
      for (j=comp2[i]; j<comp2[i+1]; j++) {
	node = comp1[j];
	for (k=node_con2[node]; k<node_con2[node+1]; k++) {
	  node2 = node_con1[k];
	  if (imap[node2] != -1) nnz++;
	}
      }
      if (nnz > max_nnz) max_nnz = nnz;
    }
    else nsub += 1;
    for (j=comp2[i]; j<comp2[i+1]; j++) imap[comp1[j]] = -1;
  }
  nsub1 = new int[nnode];
  nsub2 = new int[nsub+1]; nsub2[0] = 0;
  nsub = 0;
  int *adjnode2 = new int[max_nnode_comp+1]; adjnode2[0] = 0;
  int *adjnode1 = new int[max_nnz];
  int *part = new int[max_nnode_comp];
  int *count_parts = new int[max_nparts];
  for (i=0; i<ncomp; i++) {
    nnode_comp = comp2[i+1] - comp2[i];
    for (j=comp2[i]; j<comp2[i+1]; j++) {
      node = comp1[j];
      imap[node] = j-comp2[i];
    }
    nparts = ndof_comp[i]/ndof_target;
    if (nparts > 1) {
      nnz = 0;
      for (j=comp2[i]; j<comp2[i+1]; j++) {
	node = comp1[j];
	for (k=node_con2[node]; k<node_con2[node+1]; k++) {
	  node2 = node_con1[k];
	  if (imap[node2] != -1) {
	    adjnode1[nnz] = imap[node2];
	    nnz++;
	  }
	}
	adjnode2[j-comp2[i]+1] = nnz;
      }
      if (nparts < 8)
	metis_partgraphrecursive(&nnode_comp, adjnode2, adjnode1, vwgt,
	     adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part);
      else
	metis_partgraphkway(&nnode_comp, adjnode2, adjnode1, vwgt,
	     adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part);
      myzero(count_parts, nparts);
      for (j=0; j<nnode_comp; j++) count_parts[part[j]]++;
      for (j=0; j<nparts; j++) nsub2[nsub+j+1] = nsub2[nsub+j] + 
				     count_parts[j];
      myzero(count_parts, nparts);
      for (j=0; j<nnode_comp; j++) {
	ii = nsub2[nsub+part[j]];
	nsub1[ii+count_parts[part[j]]] = comp1[comp2[i]+j];
	count_parts[part[j]]++;
      }
      nsub += nparts;
    }
    else {
      nsub2[nsub+1] = nsub2[nsub] + nnode_comp;
      for (j=comp2[i]; j<comp2[i+1]; j++) {
	nsub1[nsub2[nsub]+j-comp2[i]] = comp1[j];
      }
      nsub += 1;
    }
    for (j=comp2[i]; j<comp2[i+1]; j++) imap[comp1[j]] = -1;
  }
  delete [] ndof_comp; delete [] adjnode1; delete [] adjnode2;
  delete [] part; delete [] count_parts;
}

void CLOP_graph::determine_overlap(int nsub1[], int nsub2[], int node_con1[], 
        int node_con2[], int nnode, int* & nosub1, int* & nosub2)
{
  int i, j, k, ii, node, node2, nnode_curr, nanode, nnz;
  memset(imap, -1, nnode*sizeof(int));
  nosub2 = new int[nsub+1]; nosub2[0] = 0;
  for (i=0; i<nsub; i++) {
    for (j=nsub2[i]; j<nsub2[i+1]; j++) {
      node = nsub1[j];
      count1[j-nsub2[i]] = node;
      imap[node] = 0;
    }
    nnode_curr = nsub2[i+1] - nsub2[i];
    for (ii=0; ii<overlap; ii++) {
      nanode = 0;
      for (j=0; j<nnode_curr; j++) {
	node = count1[j];
	if (imap[node] == 0) {
	  for (k=node_con2[node]; k<node_con2[node+1]; k++) {
	    node2 = node_con1[k];
	    if (imap[node2] == -1) {
	      count1[nnode_curr+nanode] = node2;
	      imap[node2] = 0;
	      nanode++;
	    }
	  }
	  imap[node] = 1;
	}
      }
      nnode_curr += nanode;
    }
    nosub2[i+1] = nosub2[i] + nnode_curr;
    for (j=0; j<nnode_curr; j++) imap[count1[j]] = -1;
  }
  nosub1 = new int[nosub2[nsub]];
  nnz = 0;
  for (i=0; i<nsub; i++) {
    for (j=nsub2[i]; j<nsub2[i+1]; j++) {
      node = nsub1[j];
      count1[j-nsub2[i]] = node;
      nosub1[nnz] = node;
      imap[node] = 0;
      nnz++;
    }
    nnode_curr = nsub2[i+1] - nsub2[i];
    for (ii=0; ii<overlap; ii++) {
      nanode = 0;
      for (j=0; j<nnode_curr; j++) {
	node = count1[j];
	if (imap[node] == 0) {
	  for (k=node_con2[node]; k<node_con2[node+1]; k++) {
	    node2 = node_con1[k];
	    if (imap[node2] == -1) {
	      count1[nnode_curr+nanode] = node2;
	      nosub1[nnz] = node2;
	      imap[node2] = 0;
	      nanode++; nnz++;
	    }
	  }
	  imap[node] = 1;
	}
      }
      nnode_curr += nanode;
    }
    for (j=0; j<nnode_curr; j++) imap[count1[j]] = -1;
  }
}

void CLOP_graph::determine_dof_overlap(int nosub1[], int nosub2[], int nnode,
  int dof_for_node1[], int dof_for_node2[], int* & dofpart1, int* & dofpart2)
{
  int i, j, k, node, ndof_sub, nnz;
  dofpart2 = new int[nsub+1]; dofpart2[0] = 0;
  /*
  for (i=0; i<nsub; i++) {
    cout << "MyPID = " << A->Comm().MyPID() 
	 << ", nodes for overlapping subdomain " << i << endl;
    for (j=nosub2[i]; j<nosub2[i+1]; j++) cout << nosub1[j] << " ";
    cout << endl;
  }
  */
  for (i=0; i<nsub; i++) {
    ndof_sub = 0;
    for (j=nosub2[i]; j<nosub2[i+1]; j++) {
      node = nosub1[j];
      ndof_sub += dof_for_node2[node+1] - dof_for_node2[node];
    }
    dofpart2[i+1] = dofpart2[i] + ndof_sub;
  }
  dofpart1 = new int[dofpart2[nsub]];
  nnz = 0;
  for (i=0; i<nsub; i++) {
    ndof_sub = 0;
    for (j=nosub2[i]; j<nosub2[i+1]; j++) {
      node = nosub1[j];
      for (k=dof_for_node2[node]; k<dof_for_node2[node+1]; k++) {
	dofpart1[nnz] = dof_for_node1[k];
	nnz++;
      }
    }
  }
}
