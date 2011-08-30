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
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraZoltanLib.hpp>

#ifdef HAVE_EPETRA
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#endif


#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <ctype.h>
#include <exception>

/* TODO: clean up the code */

namespace Isorropia {

#ifdef HAVE_EPETRA

namespace Epetra {

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
                               int inputType):
  Library(input_graph, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
                               Teuchos::RCP<const Epetra_MultiVector> input_coords, int inputType):
  Library(input_graph, input_coords, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			  Teuchos::RCP<CostDescriber> costs,
                          int inputType):
  Library(input_graph, costs, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
			       Teuchos::RCP<CostDescriber> costs,
			       Teuchos::RCP<const Epetra_MultiVector> input_coords,
                               Teuchos::RCP<const Epetra_MultiVector> weights, int inputType):
  Library(input_graph, costs, input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
                               int inputType):
  Library(input_matrix, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			       Teuchos::RCP<const Epetra_MultiVector> input_coords,
                               int inputType):
  Library(input_matrix, input_coords, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			  Teuchos::RCP<CostDescriber> costs, int inputType):
  Library(input_matrix, costs, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
			       Teuchos::RCP<CostDescriber> costs,
			       Teuchos::RCP<const Epetra_MultiVector> input_coords,
                               Teuchos::RCP<const Epetra_MultiVector> weights, int inputType):
  Library(input_matrix, costs, input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_MultiVector> input_coords,
                               int inputType):
  Library(input_coords, inputType)
{
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_MultiVector> input_coords,
                               Teuchos::RCP<const Epetra_MultiVector> weights,
                               int inputType):
  Library(input_coords, weights, inputType)
{
  int weightDim = weights->NumVectors();

  if (weightDim > 1){
    if (input_coords->Comm().MyPID() == 0){
      std::cout << "WARNING: Zoltan will only use the first weight of the "<< weightDim << " supplied for each object" << std::endl;
    }
  }
}

ZoltanLibClass::ZoltanLibClass(Teuchos::RCP<const Epetra_BlockMap> input_map,
                               int inputType):
  Library(input_map, inputType)
{
}



int ZoltanLibClass::precompute()
{
  std::string str1("Isorropia::ZoltanLibClass::precompute ");
  MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm default_mpicomm = MPI_COMM_WORLD;
  int itype;

  Library::precompute(); // assumes input_type_ is set

  if (input_graph_.get() || input_matrix_.get())
  {
    if (input_type_ != hgraph2d_finegrain_input_){
      computeCost();     // graph or hypergraph weights
    }
  }

  if (input_type_ == graph_input_)
    itype = ZoltanLib::QueryObject::graph_input_;
  else if (input_type_ == hgraph_input_)
    itype = ZoltanLib::QueryObject::hgraph_input_;
  else if (input_type_ == hgraph2d_finegrain_input_)
    itype = ZoltanLib::QueryObject::hgraph2d_finegrain_input_;
  else if (input_type_ == geometric_input_)
    itype = ZoltanLib::QueryObject::geometric_input_;
  else if (input_type_ == simple_input_)
    itype = ZoltanLib::QueryObject::simple_input_;
  else if (input_type_ == hgraph_graph_input_)                 // hierarchical partitioning
    itype = ZoltanLib::QueryObject::hgraph_graph_input_;
  else if (input_type_ == hgraph_geometric_input_)             // hierarchical partitioning
    itype = ZoltanLib::QueryObject::hgraph_geometric_input_;
  else if (input_type_ == graph_geometric_input_)              // hierarchical partitioning
    itype = ZoltanLib::QueryObject::graph_geometric_input_;
  else if (input_type_ == hgraph_graph_geometric_input_)       // hierarchical partitioning
    itype = ZoltanLib::QueryObject::hgraph_graph_geometric_input_;
  else
    itype = ZoltanLib::QueryObject::unspecified_input_;


  if (input_graph_.get() !=0 && input_coords_.get()!=0) //geometric and graph inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_graph_, costs_, input_coords_, weights_, itype));
#ifdef HAVE_MPI
    const  Epetra_Comm &ecomm = input_graph_->RowMap().Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
#endif
  }
  else if (input_matrix_.get() !=0 && input_coords_.get()!=0) //geometric and matrix inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_matrix_, costs_, input_coords_, weights_, itype));
#ifdef HAVE_MPI
    const Epetra_Comm &ecomm = input_matrix_->RowMatrixRowMap().Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
#endif
  }
  else if (input_graph_.get() != 0) //graph inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_graph_, costs_, itype));
#ifdef HAVE_MPI
    const  Epetra_Comm &ecomm = input_graph_->RowMap().Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
#endif
  }
  else if (input_matrix_.get() != 0) //matrix inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_matrix_, costs_, itype));
#ifdef HAVE_MPI
    const Epetra_Comm &ecomm = input_matrix_->RowMatrixRowMap().Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
#endif
  }
  else if (input_coords_.get() != 0) // coord inputs 
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_coords_, weights_));
#ifdef HAVE_MPI
    const Epetra_Comm &ecomm = input_coords_->Map().Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
#endif
  }
  else // BlockMap inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject(input_map_, itype));
#ifdef HAVE_MPI
    const  Epetra_Comm &ecomm = input_map_->Comm();
    try
    {
    const Epetra_MpiComm &empicomm = dynamic_cast<const Epetra_MpiComm &>(ecomm);
    mpicomm = empicomm.Comm();
    }
    catch (std::exception& e)
    {
        // Serial Comm with MPI
        MPI_Comm_split(default_mpicomm, ecomm.MyPID(), 0, &mpicomm);
    }
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
    if (input_type_ == graph_input_)
    {
      zoltanParamList_.set(lb_method_str, "GRAPH");
    }
    else if (input_type_ == geometric_input_)
    {
      if (!zoltanParamList_.isParameter(lb_method_str))  //MMW: Don't think this if is needed 
	zoltanParamList_.set(lb_method_str, "RCB");
    }
    else if (input_type_ == simple_input_) //not sure this is needed
    {
      zoltanParamList_.set(lb_method_str, "BLOCK");      
    }
    else if (input_type_ == hgraph_graph_input_    || input_type_ == hgraph_geometric_input_ ||
             input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_ )
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
  if (input_type_ == hgraph2d_finegrain_input_)
  {
    zoltanParamList_.set("OBJ_WEIGHT_DIM", "0");
    zoltanParamList_.set("EDGE_WEIGHT_DIM", "0");
  }
  else if (input_type_ == geometric_input_)
  {
    // We always overwrite user choice.
    // if (!zoltanParamList_.isParameter("OBJ_WEIGHT_DIM")) {
      if (weights_.get())
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "1");
      }
      else
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "0");
      }
    //}
  }
  else if(input_type_ != simple_input_) //graph or hypergraph 
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

  if (input_type_ == hgraph2d_finegrain_input_)
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

  zz_->Set_Num_Obj_Fn(ZoltanLib::QueryObject::Number_Objects, (void *)queryObject_.get());
  zz_->Set_Obj_List_Fn(ZoltanLib::QueryObject::Object_List, (void *)queryObject_.get());

  int ierr;
  num_obj_ = ZoltanLib::QueryObject::Number_Objects((void *)queryObject_.get(), &ierr);


  if (input_type_ == hgraph2d_finegrain_input_)
  {
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject::HG_CS, (void *)queryObject_.get());
  }
  if (input_type_ == hgraph_input_           || input_type_ == hgraph_graph_input_ ||
      input_type_ == hgraph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_)
  {
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject::HG_CS, (void *)queryObject_.get());
    zz_->Set_HG_Size_Edge_Wts_Fn(ZoltanLib::QueryObject::HG_Size_Edge_Weights,
				 (void *)queryObject_.get());
    zz_->Set_HG_Edge_Wts_Fn(ZoltanLib::QueryObject::HG_Edge_Weights, (void *)queryObject_.get());
  }
  if (input_type_ == graph_input_ || input_type_ == hgraph_graph_input_ ||
      input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_)
  {
    zz_->Set_Num_Edges_Multi_Fn(ZoltanLib::QueryObject::Number_Edges_Multi, (void *)queryObject_.get());
    zz_->Set_Edge_List_Multi_Fn(ZoltanLib::QueryObject::Edge_List_Multi, (void *)queryObject_.get());
  }
  if (input_type_ == geometric_input_ || input_type_ == hgraph_geometric_input_ ||
      input_type_ == graph_geometric_input_ || input_type_ == hgraph_graph_geometric_input_)
  {
    zz_->Set_Num_Geom_Fn(ZoltanLib::QueryObject::Number_Geom, (void *)queryObject_.get());
    zz_->Set_Geom_Multi_Fn(ZoltanLib::QueryObject::Geom_Multi, (void *)queryObject_.get());
  }

  return (ierr);
}

void ZoltanLibClass::computeCost()
{
  std::string str1("Isorropia::ZoltanLibClass::computeCost ");
  std::string str2;


  const Epetra_Comm &comm = input_map_->Comm();

  // If vertex/edge costs have been set, do a global operation to find
  // out how many weights were given.  Some processes may provide no
  // weights - they need to be informed that weights are being provided
  // by the application.  Do some sanity checks.

  int err = 0;
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

  if (input_graph_.get() == 0) //matrix
  {
    myNZ = input_matrix_->NumMyNonzeros();
    mySelfEdges = input_matrix_->NumMyDiagonals();
    globalNZ = input_matrix_->NumGlobalNonzeros();
    globalSelfEdges = input_matrix_->NumGlobalDiagonals();
    globalNumCols = input_matrix_->NumGlobalCols();
  }
  else //graph
  {
    myNZ = input_graph_->NumMyNonzeros();
    mySelfEdges = input_graph_->NumMyDiagonals();
    globalNZ = input_graph_->NumGlobalNonzeros();
    globalSelfEdges = input_graph_->NumGlobalDiagonals();
    globalNumCols = input_graph_->NumGlobalCols();
  }

  if (costs_.get() != 0) 
  {

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
}


void ZoltanLibClass::preCheckPartition()
{
  std::string str1("Isorropia::ZoltanLibClass::precheckPartition ");
  std::string str2;

  const Epetra_Comm &comm = input_map_->Comm();
  int localProc = comm.MyPID();
  int nprocs = comm.NumProc();

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
  int numrows = input_map_->NumGlobalElements();

  int myParts[] = {myGparts, myLparts};
  int maxParts[2];

  comm.MaxAll(myParts, maxParts, 2);

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
	// work for Epetra matrices, where rows are not owned by more than
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

int ZoltanLibClass::
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

  if (numPartSizes > 0){
    int err;
    int *wgtIdx = NULL;

    err = zz_->LB_Set_Part_Sizes(1, numPartSizes, partGIDs, wgtIdx, partSizes);

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

int ZoltanLibClass::
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

  // TODO64 - This only works if the global IDs fit in an int.

  ZOLTAN_ID_TYPE *gids=NULL;

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    gids = new ZOLTAN_ID_TYPE [num_obj_];
    if (num_obj_ && ! gids){
      throw Isorropia::Exception("Out of memory in ZoltanLibClass::color()");
      return -1;
    }

    int *intIds = queryObject_->RowMap().MyGlobalElements();
    for (int i=0; i < num_obj_; i++){
      gids[i] = (ZOLTAN_ID_TYPE)intIds[i];
    } 
  }
  else{
    gids = (ZOLTAN_ID_PTR)queryObject_->RowMap().MyGlobalElements();
  }

  int err = zz_->Color(num_gid_entries, num_obj_, gids, &properties[0]);

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    delete [] gids;
  }

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing coloring with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}

int ZoltanLibClass::
order(Teuchos::ParameterList& zoltanParamList,
      std::vector<int>& properties)
{
  zoltanParamList_ = zoltanParamList;

  precompute();

  /* Note : this works because epetra ordinal type is int */
  int num_gid_entries = 1;

  properties.resize(num_obj_);

  // TODO64 - This only works if the global IDs fit in an int.

  ZOLTAN_ID_TYPE *gids=NULL;

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    gids = new ZOLTAN_ID_TYPE [num_obj_];
    if (num_obj_ && ! gids){
      throw Isorropia::Exception("Out of memory in ZoltanLibClass::order()");
      return -1;
    }

    int *intIds = queryObject_->RowMap().MyGlobalElements();
    for (int i=0; i < num_obj_; i++){
      gids[i] = (ZOLTAN_ID_TYPE)intIds[i];
    } 
  }
  else{
    gids = (ZOLTAN_ID_PTR)queryObject_->RowMap().MyGlobalElements();
  }

  int err = zz_->Order(num_gid_entries, num_obj_, gids, &properties[0], NULL);

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    delete [] gids;
  }

  if (err != ZOLTAN_OK){
    throw Isorropia::Exception("Error computing ordering with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}



int ZoltanLibClass::postcompute()
{
  if (zz_)
    delete zz_;
  zz_ = NULL;

  return (0);
}




} // namespace EPETRA

#endif //HAVE_EPETRA

}//namespace Isorropia
