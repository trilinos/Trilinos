/*@HEADER 
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//     
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_SUPPORTGRAPH_H
#define IFPACK_SUPPORTGRAPH_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Condest.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Amesos.h"
#include "Ifpack_Condest.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_Time.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/config.hpp>

using Teuchos::RefCountPtr;
using Teuchos::rcp;
typedef std::pair<int, int> E;
using namespace boost;

typedef adjacency_list < vecS, vecS, undirectedS,
  no_property, property < edge_weight_t, double > > Graph;
typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;



template<typename T=Ifpack_Amesos> class Ifpack_SupportGraph : 
public virtual Ifpack_Preconditioner 
{

 public:
 
 //@{ \name Constructor.

 //! Constructor
 Ifpack_SupportGraph(Epetra_RowMatrix* Matrix_in);

 //@}


 //@{ \name Attribute set methods.
 //! If set true, transpose of this operator will be applied (not implemented).
 /*! This flag allows the transpose of the given operator to be used
  * implicitly.
  
   \param       
   UseTranspose_in - (In) If true, multiply by the transpose of operator,
   otherwise just use operator.
  
   \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does       
 */
 virtual int SetUseTranspose(bool UseTranspose_in);

 //@}
 

 //@{ \name Mathematical functions.
 
 //! Applies the matrix to an Epetra_MultiVector.
 /*!
    \param     
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing the result.
                          
    \return Integer error code, set to 0 if successful.  
 */
 virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

 //! Applies the preconditioner to X, returns the result in Y.
 /*!    
    \param     
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.
                            
    \return Integer error code, set to 0 if successful.    
                                                 
    \warning In order to work with AztecOO, any implementation of this method  
    must support the case where X and Y are the same object.
 */
 virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

 //! Returns the infinity norm of the global matrix (not implemented)
 virtual double NormInf() const {return(0.0);}

 //@}


 //@{ \name Attribute access functions.

 //! Returns a character string describing the operator.
 virtual const char * Label() const;

 //! Returns the current UseTranspose setting.
 virtual bool UseTranspose() const {return(UseTranspose_);};

 //! Returns true if this object can provide an approximate Inf-norm, false otherwise.
 virtual bool HasNormInf() const {return(false);};

 //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
 virtual const Epetra_Comm & Comm() const {return(Matrix_->Comm());};

 //! Returns the Epetra_Map object associated with the domain of this operator.
 virtual const Epetra_Map & OperatorDomainMap() const {return(Matrix_->OperatorDomainMap());};

 //! Returns the Epetra_Map object associated with the range of this operator.
 virtual const Epetra_Map & OperatorRangeMap() const {return(Matrix_->OperatorRangeMap());};

 //@}


 //@{ \name Construction and application methods.

 //! Returns \c true if the preconditioner has been successfully initialized.
 virtual bool IsInitialized() const
 {
   return(IsInitialized_);
 }

 //! Returns \c true if the preconditioner has been successfully computed.
 virtual bool IsComputed() const
 {
   return(IsComputed_);
 }

 //! Sets all the parameters for the preconditioner.
 /*! Parameters currently supported:
  * The input list will be copied, then passed                                       
  * to the underlying preconditioner
  *    
  * - \c "MST: forest number" : Specified the number of 
  *   times Kruskal's algorithm adds another forest to
  *   the preconditioner            
  *
  * - \c "MST: diagonal offset" : Specify the offset
  *   to add to the diagonal elements of the support
  *   graph matrix
  */
 virtual int SetParameters(Teuchos::ParameterList& List);

 //! Initialize the preconditioner    
 /*! \return 
  * 0 if successful, 1 if problems occured.
  */
 virtual int Initialize();
 //! Computes the preconditioners.
 /*! \return
  * 0 if successful, 1 if problems occurred.
  */
 virtual int Compute();

 //@}


 //@{ \name Query methods.


 //! Returns the estimated conditioner number, computes it if necessary.    
 /*!
  * not implemented    
  */
 virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
			const int MaxIters = 1550,
			const double Tol = 1e-9,
			Epetra_RowMatrix* Matrix_in = 0);

 //! Returns the computed condition number.
 virtual double Condest() const
 {
   return(Condest_);
 }

 //! Returns a const reference to the internally stored matrix.
 virtual const Epetra_RowMatrix& Matrix() const
 {
   return(*Matrix_);
 }

 //! Prints on ostream basic information about \c this object.
 virtual std::ostream& Print(std::ostream&) const;

 //! Returns the number of calls to Initialize().
 virtual int NumInitialize() const
 {
   return(NumInitialize_);
 }

 //! Returns the number of calls to Compute().
 virtual int NumCompute() const
 {
   return(NumCompute_);
 }

 //! Returns the number of calls to ApplyInverse().
 virtual int NumApplyInverse() const
 {
   return(NumApplyInverse_);
 }

 //! Returns the total time spent in Initialize().
 virtual double InitializeTime() const
 {
   return(InitializeTime_);
 }

 //! Returns the total time spent in Compute().
 virtual double ComputeTime() const
 {
   return(ComputeTime_);
 }

 //! Returns the total time spent in ApplyInverse().
 virtual double ApplyInverseTime() const
 {
   return(ApplyInverseTime_);
 }

 //! Returns the number of flops in the initialization phase.        
 virtual double InitializeFlops() const
 {
   return(InitializeFlops_);
 }

 //! Returns the total number of flops to compute the preconditioner.
 virtual double ComputeFlops() const
 {
   return(ComputeFlops_);
 }

 //! Returns the total number of flops to apply the preconditioner.
 virtual double ApplyInverseFlops() const
 {
   return(ApplyInverseFlops_);
 }


 //@}

 protected:

 // Finds the size of a subtree rooted at a vertex (work in progress)
 int treecount(const std::vector<Vertex>& v, int *subtreesize, int node);

 // Partitions the spanning tree (work in progress)
 int treepartition(int *table, int* children, int *subtreesize, std::vector<int>& roots, 
		   int node, int n, int t);

 // Finds the largest edge between two trees (work in progress)
 int largestbetween(int *table, int* children, const Graph& graph, const property_map<Graph, edge_weight_t>::type& map,
		    int tree1, int tree2, double *largest, int *extrasource, int *extratarget, int num_verts);

 // Finds a list of all the vertices of a tree rooted at input vertex (work in progress)
 int findall(std::vector<int>& tree, int root, int *table, int *children, int num_verts);


 //! Compute AMST support graph (work in progress).
 int AMST();
 
 //! Compute the support graph.
 int FindSupport();

 //! Pointers to the matrix to be preconditioned.
 Teuchos::RefCountPtr<const Epetra_RowMatrix> Matrix_;
 
 //! Pointers to the matrix of the support graph.
 Teuchos::RefCountPtr<Epetra_CrsMatrix> Support_;

 //! Contains the label of \c this object.       
 string Label_;

 //! If true, the preconditioner has been successfully initialized.     
 bool IsInitialized_;

 //! If true, the preconditioner has been successfully computed.    
 bool IsComputed_;

 //! If \c true, solve with the transpose (not supported by all solvers).     
 bool UseTranspose_;

 //! Stores a copy of the list given in SetParameters()    
 Teuchos::ParameterList List_;

 //! Contains the estimated condition number.
 double Condest_;

 //! Contains the number of successful calls to Initialize(). 
 int NumInitialize_;
 
 //! Contains the number of successful call to Compute().
 int NumCompute_;
 
 //! Contains the number of successful call to ApplyInverse().
 mutable int NumApplyInverse_;
 
 //! Contains the time for all successful calls to Initialize().
 double InitializeTime_;

 //! Contains the time for all successful calls to Compute().
 double ComputeTime_;
 
 //! Contains the time for all successful calls to ApplyInverse().
 mutable double ApplyInverseTime_;

 //! Contains the number of flops for Initialize().
 double InitializeFlops_;

 //! Contains the number of flops for Compute().
 double ComputeFlops_;

 //! Contain sthe number of flops for ApplyInverse(). 
 mutable double ApplyInverseFlops_;
 
 //! Object used for timing purposes. 
 Teuchos::RefCountPtr<Epetra_Time> Time_;

 //! Pointer to the local solver.
 Teuchos::RefCountPtr<T> Inverse_;

 //! Contains the number of forests in the support graph
 int NumForests_;

 //! Contains the Offset to add to the diagonal of the support graph
 double Offset_;

 //! Contains the option to keep the diagonal of original matrix, or weighted average
 double KeepDiag_;

}; // class Ifpack_SupportGraph<T>



//==============================================================================
template<typename T>
Ifpack_SupportGraph<T>::Ifpack_SupportGraph(Epetra_RowMatrix* Matrix_in):
Matrix_(rcp(Matrix_in,false)),
  IsInitialized_(false),
  IsComputed_(false),
  UseTranspose_(false),
  Condest_(-1.0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  InitializeFlops_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  NumForests_(1),
  Offset_(1),
  KeepDiag_(1.0)
{
  
  Teuchos::ParameterList List_in;
  SetParameters(List_in);
}
//============================================================================== 
template<typename T>
int Ifpack_SupportGraph<T>::treecount(const std::vector<Vertex>& v, int *subtreesize, int node)
{
  int parent = v[node];
 
  if(node != parent)
    {
      subtreesize[node] = subtreesize[node] + 1;
      treecount(v, subtreesize, parent);
    }

  return 0;
}
//==============================================================================
template<typename T>
int Ifpack_SupportGraph<T>::treepartition(int *table, int* children, int *subtreesize, 
					  std::vector<int>& roots,int node, int n, int t)
{
  subtreesize[node] = 1;
  for(int i = table[node]; i < table[node+1]; i++)
    {
      int child = children[i];
    
      if(subtreesize[child] > (n/t + 1))
	{
	  treepartition(table, children, subtreesize, roots, child, n, t);
	}

      if(subtreesize[child] > (n/t))
	{
	  children[i] = -children[i];
	  roots.push_back(child);
	  
	}
      else
	{
	  subtreesize[node] = subtreesize[node] + subtreesize[child];
	}
    }
  return 0;
}
//==============================================================================
template<typename T>
int Ifpack_SupportGraph<T>::largestbetween(int* table, int* children, const Graph& graph, 
					    const property_map<Graph, edge_weight_t>::type& map, 
					    int tree1, int tree2,
					   double *largest, int *extrasource, int *extratarget, int num_verts)
{
  std::vector <int> subtree1;
  std::vector <int> subtree2;

  if(tree1 < 0)
    {
      tree1 = -tree1;
    }

  if(tree2 < 0)
    {
      tree2 = -tree2;
    }

  subtree1.push_back(tree1);
  subtree2.push_back(tree2);


  findall(subtree1, tree1, table, children, num_verts);
  findall(subtree2, tree2, table, children, num_verts);


  for(int i = 0; i < subtree1.size(); i++)
    {
      for(int j = 0; j < subtree2.size(); j++)
	{

	  if(edge(subtree1[i], subtree2[j], graph).second)
	    {                                                  

	      double temp = get(map, edge(subtree1[i],subtree2[j], graph).first);                                                                                  

	      if(temp < *largest)         
		{
		  *largest = temp;
		  *extrasource = subtree1[i];
		  *extratarget = subtree2[j];
		}

	    }    
	}
    }



    return 0;
}
//==============================================================================
template<typename T>
int Ifpack_SupportGraph<T>::findall(std::vector<int>& tree, int root, int *table, int *children, int num_verts)
{
  int upper;
  if(root == num_verts-1)
    {
      upper = num_verts;
    }
  else
    {
      upper = table[root+1];
    }
  for(int i = table[root]; i < upper; i++)
    {
      int child = children[i];
     
      if(child > 0)
	{
	  tree.push_back(child);
	  findall(tree, child, table, children,num_verts);
	}
    }
}
//============================================================================== 
template<typename T>
int Ifpack_SupportGraph<T>::AMST()
{
  
 

  // Extract matrix dimensions
  int rows = (*Matrix_).NumGlobalRows();
  int cols = (*Matrix_).NumGlobalCols();
  int num_edges  = ((*Matrix_).NumMyNonzeros() - (*Matrix_).NumMyDiagonals())/2;

  // Assert square matrix
  IFPACK_CHK_ERR((rows == cols));

  // Rename for clarity
  int num_verts = rows;

  double fill = .9;
  int t = fill*pow(num_verts,.5);

  // Create data structures for the BGL code and temp data structures for extraction 
  E *edge_array = new E[num_edges];
  double *weights = new double[num_edges];
  double *shiftedweights = new double[num_edges];

  int num_entries;
  int max_num_entries = (*Matrix_).MaxNumEntries();
  double *values = new double[max_num_entries];
  int *indices = new int[max_num_entries];

  double * diagonal = new double[num_verts];
  double shift = 0;

  for(int i = 0; i < max_num_entries; i++)
    {
      values[i]=0;
      indices[i]=0;
    }

  // Extract from the epetra matrix keeping only one edge per pair (assume symmetric) 
  int k = 0;
  for(int i = 0; i < num_verts; i++)
    {
      (*Matrix_).ExtractMyRowCopy(i,max_num_entries,num_entries,values,indices);

      for(int j = 0; j < num_entries; j++)
        {

          if(i == indices[j])
            {
              diagonal[i] = values[j];
            }

          if(i < indices[j])
            {
              edge_array[k] = E(i,indices[j]);
              if(values[j] < shift)
                shift = values[j];

              weights[k] = values[j];

              k++;
            }
        }
    }

  shift = shift - 1;

  for(int i = 0; i < num_edges; i++)
    {
      shiftedweights[i] = weights[i] - shift;
    }


  std::vector<int> TreeNz(num_verts,1);

  int *TotalNz = new int[num_verts];
  for(int i = 0; i < num_verts; i++)
    {
      TotalNz[i] = 1;
    }


  Graph gtemp(edge_array, edge_array + num_edges, shiftedweights, num_verts);
  //gtemp = Graph(edge_array, edge_array + num_edges, shiftedweights, num_verts);

  property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, gtemp);

  //weightmap = get(edge_weight, gtemp);

  std::vector < Vertex > p(num_vertices(gtemp));
  prim_minimum_spanning_tree(gtemp, &p[0]);

  std::vector<int> roots;
  int numchildren[num_verts];
  int table[num_verts];
  int children[num_verts];
  int *subtreesize = new int[num_verts];
  for(std::size_t i = 0; i != p.size(); ++i)
    {
      numchildren[i] = 0;
      table[i] = 0;
      children[i] = 0;
      subtreesize[i] = 0;
    }

  for (std::size_t i = 0; i != p.size(); ++i)
    {

      if (p[i] != i)
	{
	  numchildren[p[i]] = numchildren[p[i]] + 1;
	  TreeNz[p[i]] = TreeNz[p[i]] + 1;
	  TreeNz[i] = TreeNz[i] + 1;
	  TotalNz[p[i]] = TotalNz[p[i]] + 1;
	  TotalNz[i] = TotalNz[i] + 1;
	  //std::cout << "parent[" << i << "] = " << p[i] << std::endl;
	}
      else
	{
	  roots.push_back(i);
	}
    }




  k = 0;
 
  for (std::size_t i = 0; i != p.size(); ++i)
    {
      table[i] = k;
      k = k + numchildren[i];
    }
 
 
  for (std::size_t i = 0; i != p.size(); ++i)
    {
      if(i != p[i])
	{
	  //	  std::cout << table[p[i]] << "    " << numchildren[p[i]] << "     " << table[p[i]] + numchildren[i] - 1 << std::endl;
	  children[table[p[i]] + numchildren[p[i]] - 1] = i;
	  numchildren[p[i]] = numchildren[p[i]] - 1;
	}
    }

 
  for(std::size_t i = 0; i != p.size(); ++i)
    {
      treecount(p,subtreesize, i);
    }

 
  for(int i = 0; i < roots.size(); i++)
    {
      treepartition(table, children, subtreesize, roots, roots[i], num_verts, t);
    }

  /*
  for(int i = 0; i < roots.size(); i++)
    {
      p[roots[i]] = i;
      std::cout << roots[i] << std::endl;
      }*/
  /*
  std::cout << "children" << std::endl;
  for(std::size_t i = 0; i != p.size(); ++i)
    {
      std::cout << children[i] << std::endl;
      }*/

  std::vector<int> ExtraIndices[num_verts];
  std::vector<double> ExtraValues[num_verts];

  for(int i = 0; i < num_verts; i++)
    {
      std::vector<int> temp;
      std::vector<double> temp2;
      ExtraIndices[i] = temp;
      ExtraValues[i] = temp2;
    }



  for(int i = 0; i < roots.size(); i++)
    {
      for(int j = i+1; j < roots.size(); j++)
	{
	  double largest = -shift;
	  int extrasource = -1;
	  int extratarget = -1;


	  largestbetween(table, children, gtemp, weightmap, roots[i],roots[j],&largest,&extrasource,&extratarget,num_verts);

	  if(largest < -shift)
	    {

	      if((p[extrasource] != extratarget) && (p[extratarget] != extrasource))
		{

		  ExtraIndices[extrasource].push_back(extratarget);
		  ExtraIndices[extratarget].push_back(extrasource);
		  ExtraValues[extrasource].push_back(largest+shift);
		  ExtraValues[extratarget].push_back(largest+shift);
		  TotalNz[extrasource] = TotalNz[extrasource] + 1;
		  TotalNz[extratarget] = TotalNz[extratarget] + 1;
		}
	    } 
	  
	}
	  
	  

    }


  
  Support_ = rcp(new Epetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(),TotalNz, true));


  for(int i = 0; i < num_verts; i++)
    {

      std::vector<int> Indices(TotalNz[i],0);
      std::vector<double> Values(TotalNz[i],0);
     
      int other;
      int upper;
      if(p[i] == i)
	{
	  upper = TreeNz[i] - 1;
	      
	}
      else
	{

	  std::cout << TreeNz[i] << std::endl;
	  upper = TreeNz[i] - 2;
	  Indices[TreeNz[i]-1] = p[i];

	  if(!edge(i,p[i],gtemp).second)
	    {
	      std::cout << "WTFFFFFF" << std::endl;
	    }
	  Values[TreeNz[i]-1] = get(weightmap, edge(i, p[i], gtemp).first) + shift;
	}

      for(int j = 0; j < upper; j++)
	{

	  other = children[table[i]+j];
	  if(other < 0)
	    other = -other;

	  Indices[j+1] = other;
	  Values[j+1] = get(weightmap, edge(i, other, gtemp).first) + shift;
	}

      int s = 0;
      for(int j = TreeNz[i] + 1; j < TotalNz[i]; j++)
	{
	      std::cout << "extra" << std::endl;
	      Indices[j] = ExtraIndices[i][s];
	      Values[j] = ExtraValues[i][s];

	      s++;

	}
		

      Indices[0] = i;
      Values[0] = diagonal[i];

      (*Support_).InsertGlobalValues(i,TotalNz[i],&Values[0],&Indices[0]);

    }

  (*Support_).FillComplete();

  //(*Support_).Print(std::cout);

  delete subtreesize;
  delete TotalNz;

  return 0;
}
//============================================================================== 
template<typename T>
int Ifpack_SupportGraph<T>::FindSupport()
{

  // Extract matrix dimensions                                                                  
  int rows = (*Matrix_).NumGlobalRows();
  int cols = (*Matrix_).NumGlobalCols();
  int num_edges  = ((*Matrix_).NumMyNonzeros() - (*Matrix_).NumMyDiagonals())/2;

  // Assert square matrix                                                                       
  IFPACK_CHK_ERR((rows == cols));

  // Rename for clarity                                                                         
  int num_verts = rows;

  // Create data structures for the BGL code and temp data structures for extraction            
  E *edge_array = new E[num_edges];
  double *weights = new double[num_edges];
 
  int num_entries;
  int max_num_entries = (*Matrix_).MaxNumEntries();
  double *values = new double[max_num_entries];
  int *indices = new int[max_num_entries];

  double * diagonal = new double[num_verts];

  
  for(int i = 0; i < max_num_entries; i++)
    {
      values[i]=0;
      indices[i]=0;
    }

  // Extract from the epetra matrix keeping only one edge per pair (assume symmetric)           
  int k = 0;
  for(int i = 0; i < num_verts; i++)
    {
      (*Matrix_).ExtractMyRowCopy(i,max_num_entries,num_entries,values,indices);

      for(int j = 0; j < num_entries; j++)
	{

	  if(i == indices[j])
	    {
	      diagonal[i] = values[j];
	    }

	  if(i < indices[j])
	    {
	      edge_array[k] = E(i,indices[j]);
	      weights[k] = values[j];

	      k++;
	    }
	}
    }
  
  // Create BGL graph                                                                           
  Graph g(edge_array, edge_array + num_edges, weights, num_verts);
  

  property_map < Graph, edge_weight_t >::type weight = get(edge_weight, g);

  std::vector < Edge > spanning_tree;

  // Run Kruskal, actually maximal weight ST since edges are negative                           
  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
  

  std::vector<int> NumNz(num_verts,1);

  //Find the degree of all the vertices                                                         
  for (std::vector < Edge >::iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei)
    {
      NumNz[source(*ei,g)] = NumNz[source(*ei,g)] + 1;
      NumNz[target(*ei,g)] = NumNz[target(*ei,g)] + 1;
    }
  
  
  // Create an stl vector of stl vectors to hold indices and values (neighbour edges)
  std::vector< std::vector< int > > Indices(num_verts);
  // TODO: Optimize for performance, may use arrays instead of vectors
  //std::vector<int> Indices[num_verts];
  //std::vector<double> Values[num_verts];

  std::vector< std::vector< double > > Values(num_verts);
  
  for(int i = 0; i < num_verts; i++)
    {
      std::vector<int> temp(NumNz[i],0);
      std::vector<double> temp2(NumNz[i],0);
      Indices[i] = temp;
      Values[i] = temp2;
    }
  
  int *l = new int[num_verts];
  for(int i = 0; i < num_verts; i++)
    {
      l[i] = 1;
    }
  
  // Add each spanning forest (tree) to the support graph and 
  // remove it from original graph
  for(int i = 0; i < NumForests_; i++)
    {
      if(i > 0)
	{
	  spanning_tree.clear();
	  kruskal_minimum_spanning_tree(g,std::back_inserter(spanning_tree));
	  for(std::vector < Edge >::iterator ei = spanning_tree.begin();
	      ei != spanning_tree.end(); ++ei)
	    {
	      NumNz[source(*ei,g)] = NumNz[source(*ei,g)] + 1;
	      NumNz[target(*ei,g)] = NumNz[target(*ei,g)] + 1;
	    }
	  for(int i = 0; i < num_verts; i++)
	    {
	      Indices[i].resize(NumNz[i]);
	      Values[i].resize(NumNz[i]);
	    }
	}

      for (std::vector < Edge >::iterator ei = spanning_tree.begin();
	   ei != spanning_tree.end(); ++ei)
	{
          // Assume standard Laplacian with constant row-sum.
          // Edge weights are negative, so subtract to make diagonal positive
	  Indices[source(*ei,g)][0] = source(*ei,g);
	  Values[source(*ei,g)][0] = Values[source(*ei,g)][0] - weight[*ei];
	  Indices[target(*ei,g)][0] = target(*ei,g);
	  Values[target(*ei,g)][0] = Values[target(*ei,g)][0] - weight[*ei];

	  Indices[source(*ei,g)][l[source(*ei,g)]] = target(*ei,g);
	  Values[source(*ei,g)][l[source(*ei,g)]] = weight[*ei];
	  l[source(*ei,g)] = l[source(*ei,g)] + 1;

	  Indices[target(*ei,g)][l[target(*ei,g)]] = source(*ei,g);
	  Values[target(*ei,g)][l[target(*ei,g)]] = weight[*ei];
	  l[target(*ei,g)] = l[target(*ei,g)] + 1;

	  remove_edge(*ei,g);
	}

    }

  
  // Set diagonal to weighted average of old and new values
  for(int i = 0; i < num_verts; i++)
     {
	  Values[i][0] *= (1.-KeepDiag_);
          Values[i][0] += KeepDiag_ * diagonal[i];
     }
  
  // Create the CrsMatrix for the support graph                                                 
  Support_ = rcp(new Epetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(),l, true));

 
  // Fill in the matrix with the stl vectors for each row                                       
  for(int i = 0; i < num_verts; i++)
    {
      Values[i][0] = Values[i][0] + Offset_;

      (*Support_).InsertGlobalValues(i,l[i],&Values[i][0],&Indices[i][0]);
    }
 
  (*Support_).FillComplete();

  //(*Support_).Print(std::cout);    

  delete edge_array;
  delete weights;
  delete values;
  delete indices;
  delete l;
  delete diagonal;


  return(0);
}
//==============================================================================
template<typename T>
int Ifpack_SupportGraph<T>::SetParameters(Teuchos::ParameterList& List_in)
{
  List_ = List_in;
  NumForests_ = List_in.get("MST: forest number", NumForests_);
  KeepDiag_ = List_in.get("MST: keep diagonal", KeepDiag_);
  Offset_ = List_in.get("MST: diagonal offset", Offset_); // TODO: Make this option compatible with other Ifpack preconditioners.

  return(0);
}
//==============================================================================    
template<typename T>
int Ifpack_SupportGraph<T>::Initialize()
{
  IsInitialized_ = false;
  IsComputed_ = false;

  
  if (Time_ == Teuchos::null)
    {
      Time_ = Teuchos::rcp( new Epetra_Time(Comm()) );
    }

  
  Time_->ResetStartTime();
 
  FindSupport();
  //AMST();
  Inverse_ = Teuchos::rcp(new T(Support_.get()));

  IFPACK_CHK_ERR(Inverse_->Initialize());

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();

  return(0);

}
//==============================================================================
template<typename T>
int Ifpack_SupportGraph<T>::Compute()
{
  if (IsInitialized() == false)
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();
  IsComputed_ = false;
  Condest_ = -1.0;

  IFPACK_CHK_ERR(Inverse_->Compute());

  IsComputed_ = true;                  
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();


  return(0);
}
//============================================================================== 
template<typename T>
int Ifpack_SupportGraph<T>::SetUseTranspose(bool UseTranspose_in)
{
  // store the flag -- it will be set in Initialize() if Inverse_ does not         
  // exist.   
  UseTranspose_ = UseTranspose_in;

  // If Inverse_ exists, pass it right now.                  
  if (Inverse_!=Teuchos::null)
    IFPACK_CHK_ERR(Inverse_->SetUseTranspose(UseTranspose_in));
  
  return(0);
}
//==============================================================================               
template<typename T>
int Ifpack_SupportGraph<T>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_CHK_ERR(Matrix_->Apply(X,Y));
  return(0);
}
//==============================================================================                  
template<typename T>
const char * Ifpack_SupportGraph<T>::Label() const
{
  return(Label_.c_str());
}
//==============================================================================                  
template<typename T>
int Ifpack_SupportGraph<T>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  Time_->ResetStartTime();

  Inverse_->ApplyInverse(X,Y);

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}
//==============================================================================                  
template<typename T>
std::ostream& Ifpack_SupportGraph<T>::
Print(std::ostream& os) const
{
  os << "================================================================================" << std::endl;
   os << "Ifpack_SupportGraph: " << Label () << endl << endl;
  os << "Condition number estimate = " << Condest() << endl;
  os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
  os << "Number of off diagonal entries in support graph matrix     = " << Support_->NumGlobalNonzeros()-Support_->NumGlobalDiagonals() << endl;
  os << "Fraction of off diagonals of support graph/off diagonals of original     = "
     << ((double)Support_->NumGlobalNonzeros()-Support_->NumGlobalDiagonals())/(Matrix_->NumGlobalNonzeros()-Matrix_->NumGlobalDiagonals());
  os << endl;
  os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
  os << "-----           -------   --------------       ------------     --------" << endl;
  os << "Initialize()    "   << std::setw(10) << NumInitialize_
     << "  " << std::setw(15) << InitializeTime_
     << "        0.0              0.0" << endl;
  os << "Compute()       "   << std::setw(10) << NumCompute_
     << "  " << std::setw(22) << ComputeTime_
     << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
  if (ComputeTime_ != 0.0)
    os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
  else
    os << "     " << std::setw(15) << 0.0 << endl;
  os << "ApplyInverse()  "   << std::setw(10) << NumApplyInverse_
     << "  " << std::setw(22) << ApplyInverseTime_
     << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
  if (ApplyInverseTime_ != 0.0)
    os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
  else
    os << "  " << std::setw(15) << 0.0 << endl;

  os << std::endl << std::endl;
  os << "Now calling the underlying preconditioner's print()" << std::endl;

  Inverse_->Print(os);
}
//==============================================================================
template<typename T>
double Ifpack_SupportGraph<T>::
Condest(const Ifpack_CondestType CT, const int MaxIters,
	const double Tol, Epetra_RowMatrix* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    {                
      return(-1.0);
    }
 
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);
  
  return(Condest_);
}

#endif // IFPACK_SUPPORTGRAPH_H
