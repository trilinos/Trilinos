/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//-----------------------------------------------------
// Ifpack2::ILUT is a translation of the Aztec ILUT
// implementation. The Aztec ILUT implementation was
// written by Ray Tuminaro.
// See notes below, in the Ifpack2::ILUT::Compute method.
// ABW.
//------------------------------------------------------

#ifndef IFPACK2_SUPPORTGRAPH_DEF_HPP
#define IFPACK2_SUPPORTGRAPH_DEF_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/config.hpp>
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


namespace Ifpack2 {
  


template <class MatrixType>
SupportGraph<MatrixType>::SupportGraph(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A) :
  A_ (A),
  Athresh_ (Teuchos::ScalarTraits<magnitude_type>::zero ()),
  Rthresh_ (Teuchos::ScalarTraits<magnitude_type>::one ()),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  Randomize_(1),
  NumForests_(1),
  KeepDiag_(1.0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  IsInitialized_ (false),
  IsComputed_ (false)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error,
    "Ifpack2::SupportGraph: Input matrix is null.");
}


template <class MatrixType>
SupportGraph<MatrixType>::~SupportGraph() {
}


template <class MatrixType>
void SupportGraph<MatrixType>::setParameters(const Teuchos::ParameterList& params) {
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  // Default values of the various parameters.
  magnitude_type absThresh = STM::zero ();
  magnitude_type relThresh = STM::one ();
  


  try {
    absThresh = params.get<magnitude_type> ("fact: absolute threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    absThresh = as<magnitude_type> (params.get<double> ("fact: absolute threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try {
    relThresh = params.get<magnitude_type> ("fact: relative threshold");
  }
  catch (InvalidParameterType&) {
    // Try double, for backwards compatibility.
    // The cast from double to magnitude_type must succeed.
    relThresh = as<magnitude_type> (params.get<double> ("fact: relative threshold"));
  }
  catch (InvalidParameterName&) {
    // Accept the default value.
  }

  try{
    Randomize_ = params.get<int> ("MST: randomize");
  }
  catch (InvalidParameterName&) {
  }
  
  if (absThresh != -STM::one ()) {
    Athresh_ = absThresh;
  }
  if (relThresh != -STM::one ()) {
    Rthresh_ = relThresh;
  }

  try{
    NumForests_ = params.get<int> ("MST: forest number");
  }
  catch (InvalidParameterName&) {
  }

  try{
    KeepDiag_ = params.get<double> ("MST: keep diagonal");
  }
  catch (InvalidParameterName&) {
  }

}


  // "Commit" the values only after validating all of them.  This
  // ensures that there are no side effects if this routine throws an
  // exception.

  // mfh 28 Nov 2012: The previous code would not assign Athresh_,
  // Rthresh_, RelaxValue_, or DropTolerance_ if the read-in value was
  // -1.  I don't know if keeping this behavior is correct, but I'll
  // keep it just so as not to change previous behavior.


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
SupportGraph<MatrixType>::getComm() const{
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
SupportGraph<MatrixType>::getMatrix() const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
SupportGraph<MatrixType>::getDomainMap() const
{
  return A_->getDomainMap();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
SupportGraph<MatrixType>::getRangeMap() const
{
  return A_->getRangeMap();
}


template <class MatrixType>
bool SupportGraph<MatrixType>::hasTransposeApply() const {
  return true;
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}


template <class MatrixType>
int SupportGraph<MatrixType>::getNumApply() const {
  return(NumApply_);
}


template <class MatrixType>
double SupportGraph<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}


template<class MatrixType>
double SupportGraph<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}


template<class MatrixType>
double SupportGraph<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}




template<class MatrixType>
typename SupportGraph<MatrixType>::magnitude_type
SupportGraph<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix) {

  if (! isComputed ()) {
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return(Condest_);
}


template<class MatrixType>
  void SupportGraph<MatrixType>::findSupport() {
  const scalar_type zero = STS::zero ();
  const scalar_type one = STS::one ();

  typedef std::pair<int, int> E;
  using namespace boost;
  typedef adjacency_list < vecS, vecS, undirectedS,
			 no_property, property < edge_weight_t, magnitude_type > > Graph;
  typedef typename graph_traits < Graph >::edge_descriptor Edge;
  typedef typename graph_traits < Graph >::vertex_descriptor Vertex;

  size_t rows = A_->getNodeNumRows();
  //int cols = A_->getGlobalNumCols();
  int num_edges  = (A_->getNodeNumEntries() - A_->getNodeNumDiags())/2;

  // Rename for clarity
  size_t num_verts = A_->getGlobalNumRows();

  // Create data structures for the BGL code and temp data structures for extraction
  E *edge_array = new E[num_edges];
  magnitude_type *weights = new magnitude_type[num_edges];

  size_t num_entries;
  size_t max_num_entries = A_->getNodeMaxNumRowEntries();

  scalar_type *valuestemp = new scalar_type[max_num_entries];
  

  local_ordinal_type *indicestemp = new local_ordinal_type[max_num_entries];
  
  //A_->describe(*out,Teuchos::VERB_EXTREME);

  magnitude_type * diagonal = new magnitude_type[rows];

  for(size_t i = 0; i < max_num_entries; i++)
    {
      valuestemp[i] = zero;
      indicestemp[i] = 0;
    }

  Tpetra::ArrayView<scalar_type> values(valuestemp, sizeof(scalar_type)*max_num_entries);
  Tpetra::ArrayView<local_ordinal_type> indices(indicestemp, sizeof(local_ordinal_type)*max_num_entries); 

  // Extract from the epetra matrix keeping only one edge per pair (assume symmetric) 
  int k = 0;
  for(size_t i = 0; i < rows; i++)
    {
      A_->getLocalRowCopy(i,indices,values, num_entries);
      //global_ordinal_type globalrow = A_->getRowMap()->getGlobalElement(i);
      
      for(size_t j = 0; j < num_entries; j++)
        {
	  //global_ordinal_type globalcol = A_->getRowMap()->getGlobalElement(indices[j]);
	  //std::cout << "globalrow " << globalrow << "  local row " << i << " to localcol" << indices[j] << "  to globalcol" << globalcol << " with value " << values[j] << std::endl;
	  if(i == Teuchos::as<size_t>(indices[j]))
            {
	      //std::cout << "what it's reading on the diag global entry " << globalrow << "  local entry  " << i << "  value  " << values[j] << std::endl;
              diagonal[i] = values[j];
              // Diagonal pertubation, only if requested
              if (Rthresh_)
		diagonal[i] *= Rthresh_;
              if (Athresh_)
		diagonal[i] += Athresh_;
	      //std::cout << "what diag is storing   " << diagonal[i] << std::endl;
            }

          if(i < Teuchos::as<size_t>(indices[j]))
            {
              edge_array[k] = E(i,indices[j]);
              weights[k] = values[j];
              if (Randomize_)
		{
		  // Add small random pertubation.
		  weights[k] *= (1.0 + 1e-8 * drand48());
		}

              k++;
            }
        }
    }


  // Create BGL graph
  Graph g(edge_array, edge_array + num_edges, weights, rows);

  typedef typename property_map < Graph, edge_weight_t >::type type;
  type weight = get(edge_weight, g);
  
  std::vector < Edge > spanning_tree;

  // Run Kruskal, actually maximal weight ST since edges are negative
  kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));


  
  Teuchos::ArrayRCP<size_t> NumNz(rows, 1);
  typedef typename std::vector < Edge >::iterator edge_iterator;

  //Find the degree of all the vertices
  for (edge_iterator ei = spanning_tree.begin();
       ei != spanning_tree.end(); ++ei)
    {
      LocalOrdinal localsource = source(*ei,g);
      LocalOrdinal localtarget = target(*ei,g);

      if(localsource > localtarget)
	{
	  localsource = target(*ei,g);
	  localtarget = source(*ei,g);
	}

      //NumNz[source(*ei,g)] = NumNz[source(*ei,g)] + 1;
      //NumNz[target(*ei,g)] = NumNz[target(*ei,g)] + 1;
      NumNz[localsource] = NumNz[localsource] + 1;
    }

  // Create an stl vector of stl vectors to hold indices and values (neighbour edges)
  //std::vector< std::vector< global_ordinal_type > > Indices(num_verts);
  std::vector<local_ordinal_type> **Indices = new std::vector<local_ordinal_type>*[rows];

  //std::vector< std::vector< magnitude_type > > Values(num_verts);

  std::vector<magnitude_type> **Values = new std::vector<magnitude_type>*[rows];

  for(size_t i = 0; i < rows; i++)
    {
      std::vector<local_ordinal_type> *temp = new std::vector<local_ordinal_type>(NumNz[i],0);
      std::vector<magnitude_type> *temp2 = new std::vector<scalar_type>(NumNz[i],0);
      
      Indices[i] = temp;
      Values[i] = temp2;
    }

  
  //size_t *l = new size_t[num_verts];
  Teuchos::ArrayRCP<size_t> localnumnz(rows, 1);
  for(size_t i = 0; i < rows; i++)
    {
      Indices[i]->at(0) = i;
    }
  

  // Add each spanning forest (tree) to the support graph and
  // remove it from original graph
  for(int i = 0; i < NumForests_; i++)
    {
      if(i > 0)
        {
          spanning_tree.clear();
          kruskal_minimum_spanning_tree(g,std::back_inserter(spanning_tree));
          for(edge_iterator ei = spanning_tree.begin();
              ei != spanning_tree.end(); ++ei)
            {
              NumNz[source(*ei,g)] = NumNz[source(*ei,g)] + 1;
              //NumNz[target(*ei,g)] = NumNz[target(*ei,g)] + 1;
            }
          for(size_t i = 0; i < num_verts; i++)
            {
              Indices[i]->resize(NumNz[i]);
              Values[i]->resize(NumNz[i]);
            }
        }

      for (edge_iterator ei = spanning_tree.begin();
           ei != spanning_tree.end(); ++ei)
        {
	  //LocalOrdinal localsource = A_->getRowMap()->getLocalElement(source(*ei,g));
	  //LocalOrdinal localtarget = A_->getRowMap()->getLocalElement(target(*ei,g));
	  
	  LocalOrdinal localsource = source(*ei,g);
	  LocalOrdinal localtarget = target(*ei,g);

	  if(localsource > localtarget)
	    {
	      localsource = target(*ei,g);
	      localtarget = source(*ei,g);
	    }
	  
          // Assume standard Laplacian with constant row-sum.
          // Edge weights are negative, so subtract to make diagonal positive
	  //          (*Indices[source(*ei,g)])[0] = source(*ei,g);
          Values[localtarget]->at(0) = Values[localtarget]->at(0) - weight[*ei];
          //(*Indices[target(*ei,g)])[0] = target(*ei,g);
          Values[localsource]->at(0) = Values[localsource]->at(0) - weight[*ei];

	  
          Indices[localsource]->at(localnumnz[localsource]) = localtarget;
          Values[localsource]->at(localnumnz[localsource]) = weight[*ei];
          localnumnz[localsource] = localnumnz[localsource] + 1;

          //Indices[localtarget]->at(localnumnz[localtarget]) = localsource;
          //Values[localtarget]->at(localnumnz[localtarget]) = weight[*ei];
          //localnumnz[localtarget] = localnumnz[localtarget] + 1;

          remove_edge(*ei,g);
        }

    }

  // Set diagonal to weighted average of Laplacian preconditioner
  // and the original matrix

  // First compute the "diagonal surplus" (in the original input matrix)
  // If input is a (pure, Dirichlet) graph Laplacian , this will be 0
  Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> ones(A_->getDomainMap());
  Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> surplus(A_->getRangeMap());

  ones.putScalar(one);
  A_->apply(ones, surplus);

  Teuchos::ArrayRCP<const scalar_type> surplusaccess = surplus.getData(0);

  for(size_t i = 0; i < rows; i++)
    {
      //      std::cout << "surplus  " << surplusaccess[i] << "  Values   " << Values[i][0] << "   but the diagonal is   " << diagonal[i] << std::endl;
      Values[i]->at(0) += surplusaccess[i];
      Values[i]->at(0) = KeepDiag_*diagonal[i] +
	(1.-KeepDiag_) * Values[i]->at(0);
	  //  std::cout << "new Values   " << Values[i][0] << std::endl;
    }


  // Create the CrsMatrix for the support graph
  Support_ = rcp(new Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(A_->getRowMap(),A_->getColMap(),localnumnz,Tpetra::StaticProfile));
  
  // Fill in the matrix with the stl vectors for each row
  for(size_t i = 0; i < rows; i++)
    {
      Teuchos::ArrayView<local_ordinal_type> IndicesInsert(*Indices[Teuchos::as<local_ordinal_type>(i)]);
      Teuchos::ArrayView<scalar_type> ValuesInsert(*Values[Teuchos::as<local_ordinal_type>(i)]);
      Support_->insertLocalValues(i,IndicesInsert,ValuesInsert);
    }

  Support_->fillComplete();

  
  //Support_->describe (*out, Teuchos::VERB_EXTREME);

  delete edge_array;
  delete weights;
  delete diagonal;
  delete Values;
  delete Indices;

}


template<class MatrixType>
void SupportGraph<MatrixType>::initialize() {
  Teuchos::Time timer ("SupportGraph::initialize");
  {
    Teuchos::TimeMonitor timeMon (timer);

    // clear any previous allocation
    IsInitialized_ = false;
    IsComputed_ = false;


    // check only in serial if the matrix is square
    TEUCHOS_TEST_FOR_EXCEPTION(
      getComm ()->getSize () == 1 && A_->getNodeNumRows() != A_->getNodeNumCols(),
      std::runtime_error, "Ifpack2::SuppotGraph::initialize: In serial or one-process "
      "mode, the input matrix must be square.");

    findSupport();
    
    solver = Amesos2::create<MatrixType, Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >("amesos2_cholmod", Support_);
    solver->symbolicFactorization();
  }

  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += timer.totalElapsedTime ();
}



template<class MatrixType>
void SupportGraph<MatrixType>::compute() {
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::reduceAll;

  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type> MV;


  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("SupportGraph::compute");
  { // Timer scope for timing compute()
    Teuchos::TimeMonitor timeMon (timer, true);
        
    // =================== //
    // start factorization //
    // =================== //

    
    solver->numericFactorization();


  }
  ComputeTime_ += timer.totalElapsedTime ();
  IsComputed_ = true;
  ++NumCompute_;
}


template <class MatrixType>
void SupportGraph<MatrixType>::apply(
           const Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& X,
                 Tpetra::MultiVector<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Y,
                 Teuchos::ETransp mode,
               typename MatrixType::scalar_type alpha,
               typename MatrixType::scalar_type beta) const
{
  this->template applyTempl<scalar_type,scalar_type>(X, Y, mode, alpha, beta);
}


template <class MatrixType>
template <class DomainScalar, class RangeScalar>
void SupportGraph<MatrixType>::applyTempl(
           const Tpetra::MultiVector<DomainScalar, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& X,
                 Tpetra::MultiVector<RangeScalar, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type>& Y,
                 Teuchos::ETransp mode,
               RangeScalar alpha,
               RangeScalar beta) const
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));

  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::MultiVector<DomainScalar, local_ordinal_type, global_ordinal_type, node_type> MV;

  Teuchos::Time timer ("SupportGraph::apply");
  { // Timer scope for timing apply()
    Teuchos::TimeMonitor timeMon (timer, true);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::SupportGraph::apply: You must call compute() to compute the incomplete "
      "factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::SupportGraph::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      beta != STS::zero (), std::logic_error,
      "Ifpack2::SupportGraph::apply: This method does not currently work when beta != 0.");

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      Xcopy = rcp (new MV (X));
    }
    else {
      Xcopy = rcpFromRef (X);
    }

    
    if (alpha != STS::one ()) {
      Y.scale (alpha);
    }
  
    RCP<MV> Ycopy = rcpFromRef(Y);
    
    //Xcopy->describe(*out,Teuchos::VERB_EXTREME);
    //Ycopy->describe(*out,Teuchos::VERB_EXTREME);

    solver->setB(Xcopy);
    solver->setX(Ycopy);

    solver->solve();
  }

  

  ++NumApply_;
  ApplyTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
std::string SupportGraph<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status: [initialized, computed]";
    }
    else {
      oss << "{status: [initialized, not computed]";
    }
  }
  else {
    oss << "{status: [not initialized, not computed]";
  }
  oss << ", global number of rows: " << A_->getGlobalNumRows()
      << ", global number of columns: " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}


template <class MatrixType>
void SupportGraph<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;
  Teuchos::OSTab tab (out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && getComm ()->getRank () == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << endl;
    out << "Absolute threshold = " << getAbsoluteThreshold() << endl;
    out << "Relative threshold = " << getRelativeThreshold() << endl;
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    if (isComputed()) {
      out << "Number of nonzeros in A         = " << A_->getGlobalNumEntries() << endl;
      out << "Number of edges in support graph     = " << Support_->getGlobalNumEntries()-Support_->getGlobalNumDiags() << std::endl;
      out << "Fraction of off diagonals of supportgraph/off diagonals of original      = " << (double)(Support_->getGlobalNumEntries()-Support_->getGlobalNumDiags())/((A_->getGlobalNumEntries()-A_->getGlobalNumDiags())/2) << std::endl;
      //          << " ( = " << 100.0 * (double)getGlobalNumEntries() / (double)A_->getGlobalNumEntries() << " % of A)" << endl;
      //out << "nonzeros / rows                 = " << 1.0 * getGlobalNumEntries() / U_->getGlobalNumRows() << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s) " << endl;
    out << "------------    -------    ---------------" << endl;
    out << "initialize()    " << setw(7) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << "compute()       " << setw(7) << getNumCompute()    << "    " << setw(15) << getComputeTime()    << endl;
    out << "apply()         " << setw(7) << getNumApply()      << "    " << setw(15) << getApplyTime()      << endl;
    out << "==============================================================================="                << endl;
    out << endl;

    solver->printTiming(out,verbLevel); 
  }
}


}//namespace Ifpack2

#endif /* IFPACK2_SUPPORTGRAPH_DEF_HPP */

