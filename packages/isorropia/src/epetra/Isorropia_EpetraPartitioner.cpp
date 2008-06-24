//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#include <Isorropia_EpetraPartitioner.hpp>
#ifdef HAVE_ISORROPIA_ZOLTAN
#include <Isorropia_Zoltan_Repartition.hpp>
#endif
#include <Isorropia_Exception.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#endif

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_CrsGraph> input_graph,
		  Teuchos::RefCountPtr<CostDescriber> costs,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(input_graph),
    input_matrix_(),
    costs_(costs),
    paramlist_(),
    partitioning_already_computed_(false)
{

  input_map_ = Teuchos::rcp(&(input_graph->RowMap()), false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Partitioner::
Partitioner(Teuchos::RefCountPtr<const Epetra_RowMatrix> input_matrix,
		  Teuchos::RefCountPtr<CostDescriber> costs,
                  const Teuchos::ParameterList& paramlist,
                  bool compute_partitioning_now)
  : input_map_(),
    input_graph_(),
    input_matrix_(input_matrix),
    costs_(costs),
    paramlist_(),
    partitioning_already_computed_(false)
{
  input_map_ = Teuchos::rcp(&(input_matrix->RowMatrixRowMap()),false);
  paramlist_ = paramlist;

  if (compute_partitioning_now) {
    compute_partitioning();
  }
}

Partitioner::~Partitioner()
{
}

void Partitioner::setParameters(const Teuchos::ParameterList& paramlist)
{
  paramlist_ = paramlist;
}

void Partitioner::compute_partitioning(bool force_repartitioning)
{
  int err = 0;
  std::string str1("Isorropia::Partitioner::compute_partitioning ");
  std::string str2;

  if (partitioning_already_computed_) {
    if (!force_repartitioning) {
      return;
    }
  }

  if (input_graph_.get() == 0 && input_matrix_.get() == 0) {
    str2 = "ERROR: not holding valid input graph OR matrix.";
    throw Isorropia::Exception(str1+str2);
  }

  const Epetra_Comm &comm = input_map_->Comm();
  int localProc = comm.MyPID();
  int nprocs = comm.NumProc();

  //if Isorropia was configured with Zoltan support, then we will use
  //Zoltan unless the user specified "PARTITIONING_METHOD" = "SIMPLE_LINEAR".

  bool use_zoltan = false;

  std::string partitioning_method_str("PARTITIONING_METHOD");
  std::string partitioning_method =
    paramlist_.get(partitioning_method_str, "UNSPECIFIED");

  std::string zoltan("Zoltan");

#ifdef HAVE_ISORROPIA_ZOLTAN
  if (partitioning_method != "SIMPLE_LINEAR") {
    use_zoltan = true;
  }

  if (use_zoltan){

    // Get special Zoltan parameters, if any

    Teuchos::ParameterList& sublist = paramlist_.sublist(zoltan);

    // Check that we have a valid Zoltan problem.
    //
    // Zoltan parameter LB_METHOD:
    //
    //   GRAPH:
    //     matrix must be square and symmetric
    //     if the graph has self edges (nonzeros in diagonal) we 
    //       don't flag this as an error, but we will omit them 
    //       in the Zoltan queries
    //
    //   HYPERGRAPH:
    //     matrix can be any shape, need not be symmetric
    //     this is the default if LB_METHOD is not set
     
    std::string lb_method_str("LB_METHOD");
    if (sublist.isParameter(lb_method_str)){
      std::string lb_meth = sublist.get(lb_method_str, "HYPERGRAPH");
      if (lb_meth == "GRAPH"){
        bool square = false;
        bool symmetric = false;
        if (input_graph_.get() != 0){
          if (input_graph_->NumGlobalRows() == input_graph_->NumGlobalCols()){
            square = true;
            // TODO - is there a quick way to figure out if the graph is
            // symmetric?  I can't see a way to do it.  For now we let
            // Zoltan figure this out.
            symmetric = true;
          }
        }
        else{
          if (input_matrix_->NumGlobalRows() == input_matrix_->NumGlobalCols()){
            square = true;
            // TODO - is there a quick way to figure out if the matrix is
            // symmetric?  I can't see a way to do it.  For now we let
            // Zoltan figure this out.
            symmetric = true;
          }
        } 
        if (!square){
          str2 = "LB_METHOD=GRAPH, matrix or graph must be square";
          throw Isorropia::Exception(str1+str2);
        }
        if (!symmetric){
          str2 = "LB_METHOD=GRAPH, matrix or graph must be symmetric";
          throw Isorropia::Exception(str1+str2);
        }
      }
    }
    else{
      sublist.set("LB_METHOD", "HYPERGRAPH");  // default is to do hypergraph partitioning
    }

    // Checks for Zoltan parameters NUM_GLOBAL_PARTITIONS and NUM_LOCAL_PARTITIONS.

    std::string gparts_str("NUM_GLOBAL_PARTITIONS");
    std::string lparts_str("NUM_LOCAL_PARTITIONS");
    std::string gparts("0");
    std::string lparts("0");

    if (sublist.isParameter(gparts_str)){
      gparts = sublist.get<std::string>(gparts_str);
    }
    if (sublist.isParameter(lparts_str)){
      lparts = sublist.get<std::string>(lparts_str);
    }
 
    int myGparts = atoi(gparts.c_str());
    int myLparts = atoi(lparts.c_str());
    int maxGparts, maxLparts, sumLparts;
    int fixGparts = -1;
    int fixLparts = -1;
    int numrows = input_map_->NumGlobalElements();

    comm.MaxAll(&myGparts, &maxGparts, 1);
    comm.MaxAll(&myLparts, &maxLparts, 1);

    // Fix problem if the number of rows is less than the number
    // of processes.  We need to set NUM_GLOBAL_PARTITIONS to
    // the number of rows, unless it was already set to something
    // greater than 0 but less than the number of rows.

    if (numrows < nprocs){
      if ((maxGparts == 0) || (maxGparts > numrows)){
        maxGparts = numrows;
      }
    }

    if (maxLparts > 0)
      comm.SumAll(&myLparts, &sumLparts, 1);
    else
      sumLparts = 0;

    if ((maxGparts > 0) || (maxLparts > 0)){
      // One or more processes set NGP or NLP so we need to check
      // them for validity.

      if (maxGparts != myGparts){
        fixGparts = maxGparts;    // all procs should have same NGP
      }

      if (maxGparts > 0){

        if (maxGparts > numrows){
          // This is an error because we can't split rows among partitions

          str2 = "NUM_GLOBAL_PARTITIONS exceeds number of rows (objects to be partitioned)";
          throw Isorropia::Exception(str1+str2);
        }

        if ((sumLparts > 0) && (sumLparts != maxGparts)){
          // This is an error because local partitions must sum to number of global

          str2 = "NUM_GLOBAL_PARTITIONS not equal to sum of NUM_LOCAL_PARTITIONS";
          throw Isorropia::Exception(str1+str2);
        }

        if ((sumLparts == 0) && (maxGparts < nprocs)){
          // Set NUM_LOCAL_PARTITIONS to 1 or 0, because Zoltan will divide
          // a partition across 2 or more processes when the number of
          // partitions is less than the number of processes.  This doesn't
          // work for Epetra matrices, where rows are not owned by more than
          // one process.

          fixLparts = (localProc < maxGparts) ? 1 : 0;
        }
      }
      else if (maxLparts > 0){

        // Set NUM_GLOBAL_PARTITIONS to sum of local partitions.  It's possible
        // that Zoltan does this already, but just to be safe...

        fixGparts = sumLparts;
      }

      if (fixGparts > 0){
        std::ostringstream os;
        os << fixGparts;
        std::string s = os.str();
        sublist.set(gparts_str, s);
      }
      if (fixLparts > 0){
        std::ostringstream os;
        os << fixLparts;
        std::string s = os.str();
        sublist.set(lparts_str, s);
      }
    }

    // If vertex/edge costs have been set, do a global operation to find
    // out how many weights were given.  Some processes may provide no
    // weights - they need to be informed that weights are being provided
    // by the application.  Do some sanity checks.
  
    err = 0;
    int gerr = 0;
    int base = input_map_->IndexBase();
  
    int numMyVWeights = 0;
    int numMyGWeights = 0;
    int numMyHGWeights = 0;
    int globalNumCols = 0;
    int myNZ = 0;
    int globalNZ = 0;
    int mySelfEdges = 0;
    int globalSelfEdges = 0;

    int myRows = input_map_->NumMyElements();
    int globalNumRows = input_map_->NumGlobalElements();

    if (input_graph_.get() == 0) {
      myNZ = input_matrix_->NumMyNonzeros();
      mySelfEdges = input_matrix_->NumMyDiagonals();
      globalNZ = input_matrix_->NumGlobalNonzeros();
      globalSelfEdges = input_matrix_->NumGlobalDiagonals();
      globalNumCols = input_matrix_->NumGlobalCols();
    }
    else{
      myNZ = input_graph_->NumMyNonzeros();
      mySelfEdges = input_graph_->NumMyDiagonals();
      globalNZ = input_graph_->NumGlobalNonzeros();
      globalSelfEdges = input_graph_->NumGlobalDiagonals();
      globalNumCols = input_graph_->NumGlobalCols();
    }
  
    if (costs_.get() != 0) {
  
      numMyVWeights = costs_->getNumVertices();
  
      if (costs_->haveGraphEdgeWeights()){
        for (int i=0; i<numMyVWeights; i++){
          int gid = input_map_->GID(i);
          if (gid >= base){
            numMyGWeights += costs_->getNumGraphEdges(gid);
          }
        }
      }
      numMyHGWeights = costs_->getNumHypergraphEdgeWeights();
  
      if ((numMyVWeights > 0) && (numMyVWeights != myRows)){
        str2 = "Number of my vertex weights != number of my rows";
        err = 1;
      }
      else if ((numMyGWeights > 0) && (numMyGWeights != (myNZ - mySelfEdges))){
        str2 = "Number of my graph edge weights != number of my nonzeros";
        err = 1;
      }
    }
    else{
      costs_ = Teuchos::rcp(new CostDescriber());
    }

    comm.SumAll(&err, &gerr ,1);

    if (gerr > 0){
      throw Isorropia::Exception(str1+str2);
    }
  
    int lval[4], gval[4];
    lval[0] = numMyVWeights;
    lval[1] = numMyGWeights;
    lval[2] = numMyHGWeights;
  
    comm.SumAll(lval, gval, 3);
  
    int numVWeights = gval[0];
    int numGWeights = gval[1];
    int numHGWeights = gval[2];
  
    if ((numVWeights > 0) && (numVWeights != globalNumRows)){
      str2 = "Number of vertex weights supplied by application != number of rows";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numGWeights > 0) && (numGWeights != (globalNZ - globalSelfEdges))){
      str2 = "Number of graph edge weights supplied by application != number of edges";
      throw Isorropia::Exception(str1+str2);
    }
    if ((numHGWeights > 0) && (numHGWeights < globalNumCols)){
      str2 = "Number of hyperedge weights supplied by application < number of columns";
      throw Isorropia::Exception(str1+str2);
    }

    costs_->setNumGlobalVertexWeights(numVWeights);
    costs_->setNumGlobalGraphEdgeWeights(numGWeights);
    costs_->setNumGlobalHypergraphEdgeWeights(numHGWeights);

    if (input_graph_.get() == 0) {
      err = Isorropia::Epetra::ZoltanLib::repartition(input_matrix_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }
    else {
      err = Isorropia::Epetra::ZoltanLib::repartition(input_graph_, costs_, sublist,
					  myNewElements_, exports_, imports_);
    }

    if (err != 0) {
      throw Isorropia::Exception("error in Isorropia_Zoltan::repartition");
    }
  }   // end if (use_zoltan)
#else
  if (paramlist_.isSublist(zoltan)) {
    throw Isorropia::Exception("Zoltan requested, but Zoltan not enabled.");
  }
#endif

  if (!use_zoltan){  //we'll use the built-in simple-linear partitioner.

    int nrows = input_map_->NumMyElements();

    if (nrows && costs_.get()){
      if (costs_->haveVertexWeights()){
        std::map<int, float> vwgts;
        costs_->getVertexWeights(vwgts);
  
        double *vals = new double [nrows];
  
        for (int i=0; i<nrows; i++){
          int gid = input_map_->GID(i);
          std::map<int, float>::iterator iter = vwgts.find(gid);
          if (iter == vwgts.end()){
            throw Isorropia::Exception("error 1 in simple linear repartitioning");
          }
          vals[i] = (double)iter->second;
        }
  
        weights_ = Teuchos::rcp(new Epetra_Vector(Copy, *input_map_, vals));
  
        delete [] vals;
      }
    }

    if (nrows && !weights_.get()){
      if (input_graph_.get() != 0) {
        weights_ = Teuchos::rcp(create_row_weights_nnz(*input_graph_));
      }
      else {
        weights_ = Teuchos::rcp(create_row_weights_nnz(*input_matrix_));
      }
    }

    err = Isorropia::Epetra::repartition(*input_map_,
                                               *weights_,
                                               myNewElements_,
                                               exports_, imports_);

    if (err != 0) {
      throw Isorropia::Exception("error 2 in simple linear repartitioning");
    }
  }

  partitioning_already_computed_ = true;
}

bool Partitioner::partitioning_already_computed() const
{
  return partitioning_already_computed_;
}

int Partitioner::newPartitionNumber(int myElem) const
{
  std::map<int,int>::const_iterator iter = exports_.find(myElem);
  if (iter != exports_.end()) {
    return(iter->second);
  }

  return( input_graph_->RowMap().Comm().MyPID() );
}

int Partitioner::numElemsInPartition(int partition) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("Partitioner::numElemsInPartition not implemented for non-local partitions.");
  }

  return(myNewElements_.size());
}

void
Partitioner::elemsInPartition(int partition, int* elementList, int len) const
{
  int myPart = input_map_->Comm().MyPID();
  if (partition != myPart) {
    throw Isorropia::Exception("error in Epetra_Map::MyGlobalElements");
  }

  unsigned length = len;
  if (myNewElements_.size() < length) length = myNewElements_.size();

  for(unsigned i=0; i<length; ++i) {
    elementList[i] = myNewElements_[i];
  }
}

void Partitioner::stringToUpper(std::string &s, int &changed)
{
  std::string::iterator siter;
  changed = 0;

  for (siter = s.begin(); siter != s.end() ; siter++)
  {
    if (islower(*siter)){
      *siter = toupper(*siter);
      changed++;
    }
  }
}

void Partitioner::paramsToUpper(Teuchos::ParameterList &plist, int &changed)
{
  changed = 0;

  // get a list of all parameter names in the list

  std::vector<std::string> paramNames ;
  Teuchos::ParameterList::ConstIterator pIter;

  pIter = plist.begin();

  while (1){
    //////////////////////////////////////////////////////////////////////
    // Compiler considered this while statement an error
    // while ( pIter = plist.begin() ; pIter != plist.end() ; pIter++ ){
    // }
    //////////////////////////////////////////////////////////////////////
    if (pIter == plist.end()) break;
    const std::string & nm = plist.name(pIter);
    paramNames.push_back(nm);
    pIter++;
  }

  // Change parameter names and values to upper case
  
  for (int i=0; i < paramNames.size(); i++){

    std::string origName(paramNames[i]);
    int paramNameChanged;
    stringToUpper(paramNames[i], paramNameChanged);

    if (plist.isSublist(origName)){
      Teuchos::ParameterList &sublist = plist.sublist(origName);

      int sublistChanged;
      paramsToUpper(sublist, sublistChanged);

      if (paramNameChanged){

        // this didn't work, so I need to remove the old sublist
        // and create a new one
        //
        //sublist.setName(paramNames[i]);

        Teuchos::ParameterList newlist(sublist);
        plist.remove(origName);
        plist.set(paramNames[i], newlist);
      }
    }
    else if (plist.isParameter(origName)){

      if (!plist.isType<std::string>(paramNames[i])){
        // all isorropia and zoltan parameters should be strings,
        // ignore anything else
        continue;
      }

      std::string &paramVal = plist.get<std::string>(origName);

      int paramValChanged;
      stringToUpper(paramVal, paramValChanged);

      if (paramNameChanged || paramValChanged){
        if (paramNameChanged){
          plist.remove(origName);
        }
        plist.set(paramNames[i], paramVal);
        changed++;
      }
    }
  } // next parameter or sublist
}

} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia

