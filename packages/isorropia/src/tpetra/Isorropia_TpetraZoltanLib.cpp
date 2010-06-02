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

************************************************************************
*/
//@HEADER


#ifdef HAVE_ISORROPIA_ZOLTAN

#endif
#include <Isorropia_Exception.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>


#ifdef HAVE_ISORROPIA_TPETRA
#include <Isorropia_TpetraZoltanLib.hpp>
#include <Isorropia_TpetraCostDescriber.hpp>

#ifdef HAVE_MPI

#endif //HAVE_MPI

#endif // HAVE_ISORROPIA_TPETRA


#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>

/* TODO: clean up the code */

namespace Isorropia {

#ifdef HAVE_ISORROPIA_TPETRA

namespace Tpetra {

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
                               int inputType):
  Library<Node>(input_graph, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
                               Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords, int inputType):
  Library<Node>(input_graph, input_coords, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
			  Teuchos::RCP<CostDescriber<Node> > costs,
                          int inputType):
  Library<Node>(input_graph, costs, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<int,int,Node> > input_graph,
			       Teuchos::RCP<CostDescriber<Node> > costs,
			       Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                               Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights, int inputType):
  Library<Node>(input_graph, costs, input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
                               int inputType):
  Library<Node>(input_matrix, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
			       Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                               int inputType):
  Library<Node>(input_matrix, input_coords, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
			  Teuchos::RCP<CostDescriber<Node> > costs, int inputType):
  Library<Node>(input_matrix, costs, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::RowMatrix<double,int,int,Node> > input_matrix,
			       Teuchos::RCP<CostDescriber<Node> > costs,
			       Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                               Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights, int inputType):
  Library<Node>(input_matrix, costs, input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                               int inputType):
  Library<Node>(input_coords, inputType)
{
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > input_coords,
                               Teuchos::RCP<const ::Tpetra::MultiVector<double,int,int,Node> > weights,
                               int inputType):
  Library<Node>(input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

template <typename Node>
ZoltanLibClass<Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::Map<int,int,Node> > input_map,
                               int inputType):
  Library<Node>(input_map, inputType)
{
}



template <typename Node>
int ZoltanLibClass<Node>::precompute()
{
  std::string str1("Isorropia::ZoltanLibClass::precompute ");
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  int itype;

  Library<Node>::precompute(); // assumes input_type_ is set

  if (this->input_graph_.get() || this->input_matrix_.get())
  {
    if (this->input_type_ != this->hgraph2d_finegrain_input_){
      computeCost();     // graph or hypergraph weights
    }
  }

  if (this->input_type_ == this->graph_input_)
    itype = ZoltanLib::QueryObject<Node>::graph_input_;
  else if (this->input_type_ == this->hgraph_input_)
    itype = ZoltanLib::QueryObject<Node>::hgraph_input_;
  else if (this->input_type_ == this->hgraph2d_finegrain_input_)
    itype = ZoltanLib::QueryObject<Node>::hgraph2d_finegrain_input_;
  else if (this->input_type_ == this->geometric_input_)
    itype = ZoltanLib::QueryObject<Node>::geometric_input_;
  else if (this->input_type_ == this->simple_input_)
    itype = ZoltanLib::QueryObject<Node>::simple_input_;
  else if (this->input_type_ == this->hgraph_graph_input_)                 // hierarchical partitioning
    itype = ZoltanLib::QueryObject<Node>::hgraph_graph_input_;
  else if (this->input_type_ == this->hgraph_geometric_input_)             // hierarchical partitioning
    itype = ZoltanLib::QueryObject<Node>::hgraph_geometric_input_;
  else if (this->input_type_ == this->graph_geometric_input_)              // hierarchical partitioning
    itype = ZoltanLib::QueryObject<Node>::graph_geometric_input_;
  else if (this->input_type_ == this->hgraph_graph_geometric_input_)       // hierarchical partitioning
    itype = ZoltanLib::QueryObject<Node>::hgraph_graph_geometric_input_;
  else
    itype = ZoltanLib::QueryObject<Node>::unspecified_input_;


  if (this->input_graph_.get() !=0 && this->input_coords_.get()!=0) //geometric and graph inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_graph_, this->costs_, this->input_coords_, this->weights_, itype));
#ifdef HAVE_MPI
    const  Teuchos::Comm<int> &ecomm = this->input_graph_->RowMap().Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }
  else if (this->input_matrix_.get() !=0 && this->input_coords_.get()!=0) //geometric and matrix inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_matrix_, this->costs_, this->input_coords_, this->weights_, itype));
#ifdef HAVE_MPI
    const Teuchos::Comm<int> &ecomm = input_matrix_->RowMatrixRowMap().Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }
  else if (this->input_graph_.get() != 0) //graph inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_graph_, this->costs_, itype));
#ifdef HAVE_MPI
    const  Teuchos::Comm<int> &ecomm = this->input_graph_->RowMap().Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }
  else if (this->input_matrix_.get() != 0) //matrix inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_matrix_, this->costs_, itype));
#ifdef HAVE_MPI
    const Teuchos::Comm<int> &ecomm = input_matrix_->RowMatrixRowMap().Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }
  else if (this->input_coords_.get() != 0) // coord inputs 
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_coords_, this->weights_));
#ifdef HAVE_MPI
    const Teuchos::Comm<int> &ecomm = input_coords_->Map().Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }
  else // BlockMap inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<Node>(this->input_map_, itype));
#ifdef HAVE_MPI
    const  Teuchos::Comm<int> &ecomm = input_map_->Comm();
    const Teuchos::MpiComm<int> &empicomm = dynamic_cast<const Teuchos::MpiComm<int> &>(ecomm);
    mpicomm = empicomm.Comm();
#endif
  }



  float version;
  int argcTmp=0;
  char *argvTmp[1];
  std::string lb_method_str("LB_METHOD");

  // create a Zoltan problem

  argvTmp[0] = NULL;
  Zoltan_Initialize(argcTmp, argvTmp, &version);

  zz_ = new Zoltan(mpicomm);

  if (zz_ == NULL){
    throw Isorropia::Exception("Error creating Zoltan object");
    return (-1);
  }

  //////////////////////////
  // set problem parameters
  //////////////////////////

  std::string dbg_level_str("DEBUG_LEVEL");
  if (!zoltanParamList_.isParameter(dbg_level_str)) 
  {
    zoltanParamList_.set(dbg_level_str, "0");
  }

  if (!zoltanParamList_.isParameter(lb_method_str)) //set default parameters
  {
    if (this->input_type_ == this->graph_input_)
    {
      zoltanParamList_.set(lb_method_str, "GRAPH");
    }
    else if (this->input_type_ == this->geometric_input_)
    {
      if (!zoltanParamList_.isParameter(lb_method_str))  //MMW: Don't think this if is needed 
	zoltanParamList_.set(lb_method_str, "RCB");
    }
    else if (this->input_type_ == this->simple_input_) //not sure this is needed
    {
      zoltanParamList_.set(lb_method_str, "BLOCK");      
    }
    else if (this->input_type_ == this->hgraph_graph_input_    || this->input_type_ == this->hgraph_geometric_input_ ||
             this->input_type_ == this->graph_geometric_input_ || this->input_type_ == this->hgraph_graph_geometric_input_ )
    {
      zoltanParamList_.set(lb_method_str, "HIER");
    }
    else
    {
      zoltanParamList_.set(lb_method_str, "HYPERGRAPH");
    }
  }

  // Make LB_APPROACH = PARTITION the default in Isorropia 
  std::string lb_approach_str("LB_APPROACH");
  if (!zoltanParamList_.isParameter(lb_approach_str)) {
    zoltanParamList_.set(lb_approach_str, "PARTITION");
  }

    // For fine-grain hypergraph, we don't want obj or (hyper)edge weights
  if (this->input_type_ == this->hgraph2d_finegrain_input_)
  {
    zoltanParamList_.set("OBJ_WEIGHT_DIM", "0");
    zoltanParamList_.set("EDGE_WEIGHT_DIM", "0");
  }
  else if (this->input_type_ == this->geometric_input_)
  {
    // We always overwrite user choice.
    // if (!zoltanParamList_.isParameter("OBJ_WEIGHT_DIM")) {
      if (this->weights_.get())
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "1");
      }
      else
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "0");
      }
    //}
  }
  else if(this->input_type_ != this->simple_input_) //graph or hypergraph 
  {
    if (queryObject_->haveVertexWeights()) 
    {
      if (!zoltanParamList_.isParameter("OBJ_WEIGHT_DIM")) 
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "1");
      }
    }

    if (queryObject_->haveGraphEdgeWeights() ||
        queryObject_->haveHypergraphEdgeWeights()) 
    {
      if (!zoltanParamList_.isParameter("EDGE_WEIGHT_DIM")) 
      {
        zoltanParamList_.set("EDGE_WEIGHT_DIM", "1");
      }
    }
  }

  // For fine-grain hypergraph, we will use (row, col) of nz for
  // vertex GIDs.  Don't need LIDs.

  if (this->input_type_ == this->hgraph2d_finegrain_input_)
  {
    zoltanParamList_.set("NUM_GID_ENTRIES", "2");
    zoltanParamList_.set("NUM_LID_ENTRIES", "1");
  }

  Teuchos::ParameterList::ConstIterator
    iter = zoltanParamList_.begin(),
    iter_end = zoltanParamList_.end();

  for(; iter != iter_end; ++iter) 
  {
    const std::string& name = iter->first;
    const std::string& value = Teuchos::getValue<std::string>(iter->second);
    zz_->Set_Param(name, value);
  }

  // Set the query functions

  zz_->Set_Num_Obj_Fn(ZoltanLib::QueryObject<Node>::Number_Objects, (void *)queryObject_.get());
  zz_->Set_Obj_List_Fn(ZoltanLib::QueryObject<Node>::Object_List, (void *)queryObject_.get());

  int ierr;
  num_obj_ = ZoltanLib::QueryObject<Node>::Number_Objects((void *)queryObject_.get(), &ierr);


  if (this->input_type_ == this->hgraph2d_finegrain_input_)
  {
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject<Node>::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject<Node>::HG_CS, (void *)queryObject_.get());
  }
  if (this->input_type_ == this->hgraph_input_           || this->input_type_ == this->hgraph_graph_input_ ||
      this->input_type_ == this->hgraph_geometric_input_ || this->input_type_ == this->hgraph_graph_geometric_input_)
  {
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject<Node>::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject<Node>::HG_CS, (void *)queryObject_.get());
    zz_->Set_HG_Size_Edge_Wts_Fn(ZoltanLib::QueryObject<Node>::HG_Size_Edge_Weights,
				 (void *)queryObject_.get());
    zz_->Set_HG_Edge_Wts_Fn(ZoltanLib::QueryObject<Node>::HG_Edge_Weights, (void *)queryObject_.get());
  }
  if (this->input_type_ == this->graph_input_ || this->input_type_ == this->hgraph_graph_input_ ||
      this->input_type_ == this->graph_geometric_input_ || this->input_type_ == this->hgraph_graph_geometric_input_)
  {
    zz_->Set_Num_Edges_Multi_Fn(ZoltanLib::QueryObject<Node>::Number_Edges_Multi, (void *)queryObject_.get());
    zz_->Set_Edge_List_Multi_Fn(ZoltanLib::QueryObject<Node>::Edge_List_Multi, (void *)queryObject_.get());
  }
  if (this->input_type_ == this->geometric_input_ || this->input_type_ == this->hgraph_geometric_input_ ||
      this->input_type_ == this->graph_geometric_input_ || this->input_type_ == this->hgraph_graph_geometric_input_)
  {
    zz_->Set_Num_Geom_Fn(ZoltanLib::QueryObject<Node>::Number_Geom, (void *)queryObject_.get());
    zz_->Set_Geom_Multi_Fn(ZoltanLib::QueryObject<Node>::Geom_Multi, (void *)queryObject_.get());
  }

  return (ierr);
}

template <typename Node>
void ZoltanLibClass<Node>::computeCost()
{
  std::string str1("Isorropia::ZoltanLibClass::computeCost ");
  std::string str2;


  const Teuchos::Comm<int> &comm = this->input_map_->Comm();

  // If vertex/edge costs have been set, do a global operation to find
  // out how many weights were given.  Some processes may provide no
  // weights - they need to be informed that weights are being provided
  // by the application.  Do some sanity checks.

  int err = 0;
  int gerr = 0;
  int base = this->input_map_->IndexBase();

  int numMyVWeights = 0;
  int numMyGWeights = 0;
  int numMyHGWeights = 0;
  int globalNumCols = 0;
  int myNZ = 0;
  int globalNZ = 0;
  int mySelfEdges = 0;
  int globalSelfEdges = 0;

  int myRows = this->input_map_->NumMyElements();
  int globalNumRows = this->input_map_->NumGlobalElements();

  if (this->input_graph_.get() == 0) //matrix
  {
    myNZ = this->input_matrix_->NumMyNonzeros();
    mySelfEdges = this->input_matrix_->NumMyDiagonals();
    globalNZ = this->input_matrix_->NumGlobalNonzeros();
    globalSelfEdges = this->input_matrix_->NumGlobalDiagonals();
    globalNumCols = this->input_matrix_->NumGlobalCols();
  }
  else //graph
  {
    myNZ = this->input_graph_->NumMyNonzeros();
    mySelfEdges = this->input_graph_->NumMyDiagonals();
    globalNZ = this->input_graph_->NumGlobalNonzeros();
    globalSelfEdges = this->input_graph_->NumGlobalDiagonals();
    globalNumCols = this->input_graph_->NumGlobalCols();
  }

  if (this->costs_.get() != 0) 
  {

    numMyVWeights = this->costs_->getNumVertices();

    if (this->costs_->haveGraphEdgeWeights()){
	for (int i=0; i<numMyVWeights; i++){
	int gid = this->input_map_->GID(i);
	if (gid >= base){
	  numMyGWeights += this->costs_->getNumGraphEdges(gid);
	}
	}
    }
    numMyHGWeights = this->costs_->getNumHypergraphEdgeWeights();

    if ((numMyVWeights > 0) && (numMyVWeights != myRows)){
	str2 = "Number of my vertex weights != number of my rows";
	err = 1;
    }
    else if ((numMyGWeights > 0) && (numMyGWeights != (myNZ - mySelfEdges))){
	str2 = "Number of my graph edge weights != number of my nonzeros";
	err = 1;
    }
  }
  else
  {
    this->costs_ = Teuchos::rcp(new CostDescriber<Node>());
  }

  //comm.SumAll(&err, &gerr ,1);
  reduceAll(comm, Teuchos::REDUCE_SUM, 1, &err, &gerr);


  if (gerr > 0){
    throw Isorropia::Exception(str1+str2);
  }

  int lval[4], gval[4];
  lval[0] = numMyVWeights;
  lval[1] = numMyGWeights;
  lval[2] = numMyHGWeights;

  //comm.SumAll(lval, gval, 3);
  reduceAll(comm, Teuchos::REDUCE_SUM, 3, lval, gval);

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

  this->costs_->setNumGlobalVertexWeights(numVWeights);
  this->costs_->setNumGlobalGraphEdgeWeights(numGWeights);
  this->costs_->setNumGlobalHypergraphEdgeWeights(numHGWeights);
}

template <typename Node>
void ZoltanLibClass<Node>::preCheckPartition()
{
  std::string str1("Isorropia::ZoltanLibClass::precheckPartition ");
  std::string str2;

  const Teuchos::Comm<int> &comm = this->input_map_->Comm();
  int localProc = Teuchos::rank<int>(comm);
  int nprocs = Teuchos::size<int>(comm);

  // Checks for Zoltan parameters NUM_GLOBAL_PARTS and NUM_LOCAL_PARTS.

  std::string gparts_str("NUM_GLOBAL_PARTS");
  std::string lparts_str("NUM_LOCAL_PARTS");
  std::string gparts("0");
  std::string lparts("0");

  if (zoltanParamList_.isParameter(gparts_str)){
    gparts = zoltanParamList_.get<std::string>(gparts_str);
  }
  if (zoltanParamList_.isParameter(lparts_str)){
    lparts = zoltanParamList_.get<std::string>(lparts_str);
  }

  int myGparts = atoi(gparts.c_str());
  int myLparts = atoi(lparts.c_str());
  int maxGparts, maxLparts, sumLparts;
  int fixGparts = -1;
  int fixLparts = -1;
  int numrows = this->input_map_->NumGlobalElements();

  int myParts[] = {myGparts, myLparts};
  int maxParts[2];

  //comm.MaxAll(myParts, maxParts, 2);
  reduceAll(comm, Teuchos::REDUCE_SUM, 2, myParts, maxParts);


  maxGparts = maxParts[0];
  maxLparts = maxParts[1];

  // Fix problem if the number of rows is less than the number
  // of processes.  We need to set NUM_GLOBAL_PARTS to
  // the number of rows, unless it was already set to something
  // greater than 0 but less than the number of rows.

  if (numrows < nprocs){
    if ((maxGparts == 0) || (maxGparts > numrows)){
      maxGparts = numrows;
    }
  }

  if (maxLparts > 0)
  {
    //comm.SumAll(&myLparts, &sumLparts, 1);
    reduceAll(comm, Teuchos::REDUCE_SUM, 1, &myLparts, &sumLparts);
  }
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
	str2 = "NUM_GLOBAL_PARTS exceeds number of rows (objects to be partitioned)";
	throw Isorropia::Exception(str1+str2);
      }

      if ((sumLparts > 0) && (sumLparts != maxGparts)){
	// This is an error because local partitions must sum to number of global
	str2 = "NUM_GLOBAL_PARTS not equal to sum of NUM_LOCAL_PARTS";
	throw Isorropia::Exception(str1+str2);
      }

      if ((sumLparts == 0) && (maxGparts < nprocs)){
	// Set NUM_LOCAL_PARTS to 1 or 0, because Zoltan will divide
	// a partition across 2 or more processes when the number of
	// partitions is less than the number of processes.  This doesn't
	// work for Tpetra matrices, where rows are not owned by more than
	// one process.

	fixLparts = (localProc < maxGparts) ? 1 : 0;
      }
    }
    else if (maxLparts > 0){

      // Set NUM_GLOBAL_PARTS to sum of local partitions.  
      // Zoltan does this already, but just to be safe...

      fixGparts = sumLparts;
    }
    if (fixGparts > 0){
      std::ostringstream os;
      os << fixGparts;
      std::string s = os.str();
      zoltanParamList_.set(gparts_str, s);
    }
    if (fixLparts >= 0){
      std::ostringstream os;
      os << fixLparts;
      std::string s = os.str();
      zoltanParamList_.set(lparts_str, s);
    }
  }
}

template <typename Node>
int ZoltanLibClass<Node>::
repartition(Teuchos::ParameterList& zoltanParamList,
	    std::vector<int>& properties,
	    int& exportsSize,
	    std::vector<int>& imports)
{
  zoltanParamList_  = zoltanParamList;

  // Avoid to construct import list.
  // Perhaps "PARTITION ASSIGNMENTS" will be better in term of performance.
  zoltanParamList_.set("RETURN_LISTS", "EXPORT AND IMPORT");

  preCheckPartition();
  precompute();

  // Set part sizes

  if (this->numPartSizes > 0)
  {
    int err;
    int *wgtIdx = NULL;

    err = zz_->LB_Set_Part_Sizes(1, this->numPartSizes, this->partGIDs, wgtIdx, this->partSizes);

    if (err != ZOLTAN_OK){
      throw Isorropia::Exception("Error in LB_Set_Part_Sizes");
      return -1;
    }
  }

  //Generate Load Balance
  int changes=0, num_gid_entries=0, num_lid_entries=0, num_import=0, num_export=0;
  ZOLTAN_ID_PTR import_global_ids=NULL, import_local_ids=NULL;
  ZOLTAN_ID_PTR export_global_ids=NULL, export_local_ids=NULL;
  int * import_procs=NULL, * export_procs=NULL;
  int *import_to_part=NULL, *export_to_part=NULL;

  int err = zz_->LB_Partition(changes, num_gid_entries, num_lid_entries,
   num_import, import_global_ids, import_local_ids, import_procs, import_to_part,
   num_export, export_global_ids, export_local_ids, export_procs, export_to_part );

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing partitioning with Zoltan");
    return -1;
  }

  exportsSize = num_export;
  imports.clear();
  imports.assign(import_global_ids, import_global_ids + num_import);

  properties.assign(num_obj_, queryObject_->RowMap().Comm().MyPID());

  for( int i = 0; i < num_export; ++i ) 
  {
    properties[export_local_ids[i]] = export_to_part[i];
  }

  //Free Zoltan Data
  zz_->LB_Free_Part(&import_global_ids, &import_local_ids,
		     &import_procs, &import_to_part);
  zz_->LB_Free_Part(&export_global_ids, &export_local_ids,
		     &export_procs, &export_to_part);

  postcompute();

  return (0);
}

template <typename Node>
int ZoltanLibClass<Node>::
color(Teuchos::ParameterList& zoltanParamList,
      std::vector<int>& properties)
{
  zoltanParamList_ = zoltanParamList;
  // Distance 2 coloring by default.
  if (!zoltanParamList_.isParameter("COLORING_PROBLEM"))
    zoltanParamList_.set("COLORING_PROBLEM", "DISTANCE-2");

  precompute();


  //Generate Load Balance
  int  num_gid_entries = 1;


  // Use ColMap to have directly answers for columns !
  properties.resize(num_obj_);
  int err = zz_->Color(num_gid_entries, num_obj_,
 		      (ZOLTAN_ID_PTR)queryObject_->RowMap().MyGlobalElements(), &properties[0]);

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing coloring with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}

template <typename Node>
int ZoltanLibClass<Node>::
order(Teuchos::ParameterList& zoltanParamList,
      std::vector<int>& properties)
{
  zoltanParamList_ = zoltanParamList;

  precompute();

  /* Note : this works because tpetra<int,int> ordinal type is int */
  int num_gid_entries = 1;

  properties.resize(num_obj_);
  int err = zz_->Order(num_gid_entries, num_obj_,
		       (ZOLTAN_ID_PTR)queryObject_->RowMap().MyGlobalElements(), &properties[0], NULL);

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing ordering with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}


template <typename Node>
int ZoltanLibClass<Node>::postcompute()
{
  if (zz_)
    delete zz_;
  zz_ = NULL;

  return (0);
}




} // namespace TPETRA

#endif //HAVE_ISORROPIA_TPETRA

}//namespace Isorropia
