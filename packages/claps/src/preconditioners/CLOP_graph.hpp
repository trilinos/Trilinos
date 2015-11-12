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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef CLOP_GRAPH_HPP
#define CLOP_GRAPH_HPP
#include "Epetra_ConfigDefs.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_IntVector.h"

extern "C"{
  void metis_partgraphrecursive(int* n, int xadj[], int adjncy[], 
      int vwgt[], int adjwgt[], int* wgtflag, int* numflag, int* nparts,
      int options[], int* edgecut, int part[]);
  void metis_partgraphkway(int* n, int xadj[], int adjncy[], 
      int vwgt[], int adjwgt[], int* wgtflag, int* numflag, int* nparts,
      int options[], int* edgecut, int part[]);
}

class CLOP_graph 
{
 public: // functions
  CLOP_graph(
     const Epetra_CrsMatrix* A_,   // stiffness matrix
     const Epetra_IntVector* ND_,  // nodes for dofs
     const int overlap_,           // overlap
     const int partition_option_,  // partitioning option
     const int atype_,             // analysis type
     const int ndim_,              // spatial dimension
     int iwork1[],                 // integer work array
     int iwork2[],                 // integer work array
     int* & dofpart1,              // dofpart1[dofpart2[i]:dofpart2[i+1]-1] =
     int* & dofpart2,              //  dofs in overlapping subdomain i
     int & npart);                 // number of subdomains for processor
  ~CLOP_graph();
 private: // functions
  void construct_node_graph(int* & node_con1, int* & node_con2,
	     int* & dof_for_node1, int* & dof_for_node2, int & nnode);
  void determine_components(int A1[], int A2[], int N, 
			    int* & comp1, int* & comp2, int & ncomp);
  void determine_subdomains(int comp1[], int comp2[], int ncomp, 
	     int node_con1[], int node_con2[], int dof_for_node2[], 
	     int nnode, int* & nsub1, int* & nsub2);
  void determine_overlap(int nsub1[], int nsub2[], int node_con1[], 
	     int node_con2[], int nnode, int* &  nosub1, int* & nosub2);
  void determine_dof_overlap(int nosub1[], int nosub2[], int nnode, 
       int dof_for_node1[], int dof_for_node2[], int* & dofpart1, 
       int* & dofpart2);

 private: // variables
  const Epetra_CrsMatrix *A;
  const Epetra_IntVector *ND;
  const int overlap, partition_option, atype, ndim;
  int *count1, *imap;

  int ndof, nsub, ndof_target;
};
#endif // CLOP_GRAPH_HPP
