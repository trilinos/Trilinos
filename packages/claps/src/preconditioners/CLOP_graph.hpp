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
