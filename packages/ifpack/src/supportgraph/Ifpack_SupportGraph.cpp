#include "Ifpack_Preconditioner.h"
#include "Ifpack_SupportGraph.h"
#include "Ifpack_Condest.h"
#include "Ifpack_ConfigDefs.h"

#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Util.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"


#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>



using Teuchos::RefCountPtr;
using Teuchos::rcp;
typedef std::pair<int, int> E;
using namespace boost;

typedef adjacency_list < vecS, vecS, undirectedS,
			 no_property, property < edge_weight_t, int > > Graph;
typedef graph_traits < Graph >::edge_descriptor Edge;
typedef graph_traits < Graph >::vertex_descriptor Vertex;


//==============================================================================
Ifpack_SupportGraph::Ifpack_SupportGraph(Epetra_RowMatrix* A):
  A_(rcp(A,false)),
  Comm_(A->Comm()),
  UseTranspose_(false),
  Condest_(-1.0),
  NumForests_(1),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0),
  ApplyInverseTime_(0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0)
{
  Teuchos::ParameterList List;
  SetParameters(List);
  Problem_ = Teuchos::rcp( new Epetra_LinearProblem );
}
//==============================================================================
Ifpack_SupportGraph::~Ifpack_SupportGraph()
{
}
//==============================================================================
int Ifpack_SupportGraph::SetParameters(Teuchos::ParameterList& List_in)
{
  List_ = List_in;
  Label_ = List_in.get("amesos: solver type", Label_);
  NumForests_ = List_in.get("MST: forest number", NumForests_);
  return(0);
}
//==============================================================================
int Ifpack_SupportGraph::Initialize()
{
  IsInitialized_ = false;
  if(Time_ == Teuchos::null)
    Time_ = Teuchos::rcp(new Epetra_Time(Comm()));
    
  // Check #procs, won't work in parallel right now
  if(!Comm().MyPID())
    {
      // Extract matrix dimensions
      int rows = (*A_).NumGlobalRows();
      int cols = (*A_).NumGlobalCols();
      int num_edges  = ((*A_).NumGlobalNonzeros() - (*A_).NumGlobalDiagonals())/2;
      
     
      // Assert square matrix
      IFPACK_CHK_ERR((rows == cols));

      // Rename for clarity
      int num_verts = rows;

      // Create data structures for the BGL code and temp data structures for extraction
      E *edge_array = new E[num_edges];
      double *weights = new double[num_edges];

      int num_entries;
      double *values = new double[num_verts];
      int *indices = new int[num_verts];
      for(int i = 0; i < num_verts; i++)
	{
	  values[i]=0;
	  indices[i]=0;
	}
      
      // Extract from the epetra matrix keeping only one edge per pair (assume symmetric)
      int k = 0;
      for(int i = 0; i < num_verts; i++)
	{

	  (*A_).ExtractMyRowCopy(i,num_verts,num_entries,values,indices);

	  for(int j = 0; j < num_entries; j++)
	    {
	      if(i > indices[j])
		{

		  edge_array[k] = E(i+1,indices[j]+1);
		  double temp = values[j];
		  weights[k] = temp;

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
	  NumNz[source(*ei,g)-1] = NumNz[source(*ei,g)-1] + 1;
	  NumNz[target(*ei,g)-1] = NumNz[target(*ei,g)-1] + 1;
	}


      // Create an array of stl vectors to hold indices and values (neighbour edges)

      std::vector<int> Indices[num_verts];
      std::vector<double> Values[num_verts];

      
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

      for(int i = 0; i < NumForests_; i++)
	{
	  
	  if(i > 0)
	    {
	      spanning_tree.clear();
	      kruskal_minimum_spanning_tree(g,std::back_inserter(spanning_tree));
	      for(std::vector < Edge >::iterator ei = spanning_tree.begin();
		  ei != spanning_tree.end(); ++ei)
		{
		  NumNz[source(*ei,g)-1] = NumNz[source(*ei,g)-1] + 1;
		  NumNz[target(*ei,g)-1] = NumNz[target(*ei,g)-1] + 1;
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
	 
	      Values[source(*ei,g)-1][0] = Values[source(*ei,g)-1][0] - weight[*ei];
	      Values[target(*ei,g)-1][0] = Values[target(*ei,g)-1][0] - weight[*ei];
	      Indices[source(*ei,g)-1][0] = source(*ei,g)-1;
	      Indices[target(*ei,g)-1][0] = target(*ei,g)-1;

	      Indices[source(*ei,g)-1][l[source(*ei,g)-1]] = target(*ei,g)-1;
	      Values[source(*ei,g)-1][l[source(*ei,g)-1]] = weight[*ei];
	      l[source(*ei,g)-1] = l[source(*ei,g)-1] + 1;

	      Indices[target(*ei,g)-1][l[target(*ei,g)-1]] = source(*ei,g)-1;
	      Values[target(*ei,g)-1][l[target(*ei,g)-1]] = weight[*ei];
	      l[target(*ei,g)-1] = l[target(*ei,g)-1] + 1;

	      remove_edge(*ei,g);
	    }
	}

      // Create the CrsMatrix for the support graph
      B_ = rcp(new Epetra_CrsMatrix(Copy, Matrix().RowMatrixRowMap(),l, true));
      

      // Fill in the matrix with the stl vectors for each row
      for(int i = 0; i < num_verts; i++)
	{
	  Values[i][0] = Values[i][0] + 1;
	  (*B_).InsertGlobalValues(i,l[i],&Values[i][0],&Indices[i][0]);
	}
      
      (*B_).FillComplete();
      
      //(*B_).Print(std::cout);


      
      delete edge_array;
      delete weights;
      delete values;
      delete indices;
      delete l;

      // Create the Amesos to factor the system
      Problem_->SetOperator(const_cast<Epetra_CrsMatrix*>(B_.get()));
      
      Amesos Factory;
      Solver_ = Teuchos::rcp( Factory.Create("Amesos_Klu",*Problem_) );

      if(Solver_ == Teuchos::null)
	{
	  //try to create KLU, it is generally enabled
	  Solver_ = Teuchos::rcp(Factory.Create("Amesos_Klu",*Problem_));
	}

      if (Solver_ == Teuchos::null)
	IFPACK_CHK_ERR(-1);

      IFPACK_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_));
      Solver_->SetParameters(List_);
      IFPACK_CHK_ERR(Solver_->SymbolicFactorization());
      
    }
  else
    {
      std::cout << "Not yet setup for multiple processors" << std::endl;
    }


  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();

  return(0);
}
//==============================================================================
int Ifpack_SupportGraph::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  IsComputed_ = false;
  Time_->ResetStartTime();

  IFPACK_CHK_ERR(Solver_->NumericFactorization());
  
  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();

  return(0);
}
//============================================================================== 
int Ifpack_SupportGraph::Apply(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const
{
  return(0);
}
//==============================================================================
int Ifpack_SupportGraph::ApplyInverse(const Epetra_MultiVector& X,
				      Epetra_MultiVector& Y) const
{
  if (IsComputed() == false)
    IFPACK_CHK_ERR(-1);


  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-1); // wrong input 

  Time_->ResetStartTime();


  Teuchos::RefCountPtr<const Epetra_MultiVector> Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = Teuchos::rcp( new Epetra_MultiVector(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  Problem_->SetLHS(&Y);
  Problem_->SetRHS((Epetra_MultiVector*)Xcopy.get());
  IFPACK_CHK_ERR(Solver_->Solve());

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();

  return(0);
}
//==============================================================================
const char* Ifpack_SupportGraph::Label() const
{
  return((char*)Label_.c_str());
}
//==============================================================================              
int Ifpack_SupportGraph::SetUseTranspose(bool UseTranspose_in)
{
  // store the value in UseTranspose_. If we have the solver, we pass to it                   
  // right away, otherwise we wait till when it is created.                                   
  UseTranspose_ = UseTranspose_in;
  if (Solver_ != Teuchos::null)
    IFPACK_CHK_ERR(Solver_->SetUseTranspose(UseTranspose_in));

  return(0);
}
//==============================================================================
double Ifpack_SupportGraph::Condest(const Ifpack_CondestType CT,
				    const int MaxIters, const double Tol,
				    Epetra_RowMatrix* Matrix_in)
{
  if(!IsComputed())
    return(-1.0);

  if(Condest_ == -1.0)
    Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}
//==============================================================================
std::ostream&
Ifpack_SupportGraph::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_SupportGraph: " << Label () << endl << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_
       << "  " << std::setw(15) << InitializeTime_
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
      os << endl;
  }


 return(os);
}
//============================================================================== 
