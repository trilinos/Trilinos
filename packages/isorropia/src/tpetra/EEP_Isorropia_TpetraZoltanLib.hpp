//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
//@HEADER

#ifndef _Isorropia_TpetraZoltanLib_hpp_
#define _Isorropia_TpetraZoltanLib_hpp_

//#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <EEP_Isorropia_TpetraLibrary.hpp>
#include <EEP_Isorropia_TpetraCostDescriber.hpp>

#include <EEP_QueryObject.hpp>
#include <zoltan_cpp.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_CrsGraph_decl.hpp>
#include <exception>

#ifdef HAVE_MPI
#include <Teuchos_DefaultComm.hpp>
#endif

namespace Isorropia {

namespace Tpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class ZoltanLibClass : public Isorropia::Tpetra::Library<LocalOrdinal, GlobalOrdinal, Node> {
public:
  ZoltanLibClass(Teuchos::RCP< const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > & input_graph, // EEP__
                 Teuchos::RCP< CostDescriber<LocalOrdinal, GlobalOrdinal, Node> > & costs, int inputType=ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::unspecified_input_);

  /** Method to partition the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern partitioning. 

      \param[out]  newPartitions The new partition for each of my objects, in
                                   local ID order.  The objects may be rows or
                               non-zeroes (for
                               CrsGraph and RowMatrix input) or coordinates (for
                               MultiVector input).  Partition numbers can range from
                               zero to numProcs-1.
      \param[out]  exportsSize  The number of my objects that will be exported to
                              another process under the new partitioning.  This is
                             also the number of elements in newPartitions that are
                             not equal to my process rank.
      \param[out]  imports   A list of the global IDs of the objects that will be
                            imported to my process under the new partitioning
   */

  virtual int
  repartition(Teuchos::ParameterList& paramlist,
              std::vector<int>& newPartitions,
              int& exportsSize,
              std::vector<int>& imports);

  /** Method to color the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern coloring. 

      \param[out]  colorAssignment A list of integers indicating the coloring of
                              the object, in local ID order.
  */
  virtual int
  color(Teuchos::ParameterList& paramlist,
        std::vector<int>& colorAssignment);

  /** Method to order the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern ordering . 

      \param[out]  orderAssignment A list of integers indicating the ordering of
                              the object, in local ID order.
  */
  virtual int
  order(Teuchos::ParameterList& paramlist,
        std::vector<int>& orderAssignment);

protected:
  virtual int precompute();
  virtual int postcompute();
  void computeCost();
  void preCheckPartition();

  void setParameterList(Teuchos::ParameterList& zoltanParamList);

private:
  Teuchos::ParameterList zoltanParamList_;
  std::string partMethod_; // stores partitioning method used, perhaps should be in EpetraLibrary?
  Zoltan *zz_;
  Teuchos::RCP< ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node> > queryObject_;
  int num_obj_;

};//class ZoltanLibClass

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::ZoltanLibClass(Teuchos::RCP<const ::Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>> & input_graph, // EEP__
                          Teuchos::RCP<CostDescriber<LocalOrdinal, GlobalOrdinal, Node>> & costs,
                          int inputType):
  Library<LocalOrdinal, GlobalOrdinal, Node>(input_graph, costs, inputType)
{
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::precompute()
{
  std::cout << "EEP Entering isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::precompute()"
            << ": costs_ = " << this->costs_
	    << std::endl;

  std::string str1("Isorropia::ZoltanLibClass::precompute ");
  MPI_Comm mpicomm = zoltan_get_global_comm();
#ifdef HAVE_MPI
  //MPI_Comm default_mpicomm = zoltan_get_global_comm(); // EEP__ unused by compiler
#endif // HAVE_MPI
  int itype;

  Library<LocalOrdinal, GlobalOrdinal, Node>::precompute(); // assumes input_type_ is set

  if (this->input_graph_.get()) //  || input_matrix_.get()) // EEP
  {
    if (this->input_type_ != this->hgraph2d_finegrain_input_){
      computeCost();     // graph or hypergraph weights
    }
  }

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::precompute(), pos 000"
            << ": costs_ = " << this->costs_
	    << std::endl;

  if (this->input_type_ == this->graph_input_)
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::graph_input_;
  else if (this->input_type_ == this->hgraph_input_)
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::hgraph_input_;
  else if (this->input_type_ == this->hgraph2d_finegrain_input_)
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::hgraph2d_finegrain_input_;
  else if (this->input_type_ == this->geometric_input_)
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::geometric_input_;
  else if (this->input_type_ == this->simple_input_)
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::simple_input_;
  else if (this->input_type_ == this->hgraph_graph_input_)                 // hierarchical partitioning
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::hgraph_graph_input_;
  else if (this->input_type_ == this->hgraph_geometric_input_)             // hierarchical partitioning
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::hgraph_geometric_input_;
  else if (this->input_type_ == this->graph_geometric_input_)              // hierarchical partitioning
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::graph_geometric_input_;
  else if (this->input_type_ == this->hgraph_graph_geometric_input_)       // hierarchical partitioning
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::hgraph_graph_geometric_input_;
  else
    itype = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::unspecified_input_;

  if (this->input_graph_.get() != 0) //graph inputs
  {
    queryObject_ =  Teuchos::rcp(new ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>(this->input_graph_, this->costs_, itype));
#ifdef HAVE_MPI
#if 0 // EEP___
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
#endif
  }
  else
  {
    throw std::runtime_error("Invalid situation in ZoltanLibClass<>::precompute()");
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
    throw std::runtime_error/*Isorropia::Exception*/("Error creating Zoltan object");
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

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::precompute(), pos 001"
            << ": this->input_type_ = " << this->input_type_
    //<< ", queryObject_->haveVertexWeights() = " << queryObject_->haveVertexWeights()
    //<< ", queryObject_->haveGraphEdgeWeights() = " << queryObject_->haveGraphEdgeWeights()
    //<< ", queryObject_->haveHypergraphEdgeWeights() = " << queryObject_->haveHypergraphEdgeWeights()
	    << std::endl;
#if 0 // EEP
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
      if (false) // this->weights_.get()) // EEP
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "1");
      }
      else
      {
        zoltanParamList_.set("OBJ_WEIGHT_DIM", "0");
      }
    //}
  }
  else if(this->input_type_ != this->simple_input_) //graph or hypergraph // AquiAquiAqui
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
#endif // EEP
  
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

  zz_->Set_Num_Obj_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::Number_Objects, (void *)queryObject_.get());
  zz_->Set_Obj_List_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::Object_List, (void *)queryObject_.get());

  int ierr;
  num_obj_ = ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::Number_Objects((void *)queryObject_.get(), &ierr);
  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::precompute(), pos 002"
            << ": num_obj_ = " << num_obj_
	    << std::endl;

  if (this->input_type_ == this->hgraph_input_           || this->input_type_ == this->hgraph_graph_input_ ||
      this->input_type_ == this->hgraph_geometric_input_ || this->input_type_ == this->hgraph_graph_geometric_input_)
  {
    zz_->Set_HG_Size_CS_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Size_CS, (void *)queryObject_.get());
    zz_->Set_HG_CS_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_CS, (void *)queryObject_.get());
    zz_->Set_HG_Size_Edge_Wts_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Size_Edge_Weights,
                                 (void *)queryObject_.get());
    zz_->Set_HG_Edge_Wts_Fn(ZoltanLib::QueryObject<LocalOrdinal, GlobalOrdinal, Node>::HG_Edge_Weights, (void *)queryObject_.get());
  }

  std::cout << "EEP Leaving isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::precompute()" << std::endl;
  return (ierr);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void
ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::ZoltanLibClass::computeCost()
{
  std::cout << "EEP Entering isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::computeCost()" << std::endl;
  std::string str1("Isorropia::ZoltanLibClass::computeCost ");
  std::string str2;


  //const Epetra_Comm &comm = input_map_->Comm(); // EEP

  // If vertex/edge costs have been set, do a global operation to find
  // out how many weights were given.  Some processes may provide no
  // weights - they need to be informed that weights are being provided
  // by the application.  Do some sanity checks.

  //int err = 0; // EEP
  int gerr = 0;
  //int base = this->input_map_->getIndexBase(); // EEP

  int numMyVWeights = 0;
  int numMyGWeights = 0;
  int numMyHGWeights = 0;
  int globalNumCols = 0;
  int myNZ = 0;
  int globalNZ = 0;
  int mySelfEdges = 0;
  int globalSelfEdges = 0;

  //int myRows = this->input_map_->getLocalNumElements(); // EEP
  int globalNumRows = this->input_map_->getGlobalNumElements();

  if (this->input_graph_.get() == 0) //matrix
  {
    throw std::runtime_error("EEP null input_graph_ pointer in isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::computeCost()");
#if 0 // EEP
    myNZ = input_matrix_->NumMyNonzeros();
    mySelfEdges = input_matrix_->NumMyDiagonals();
    globalNZ = input_matrix_->NumGlobalNonzeros();
    globalSelfEdges = input_matrix_->NumGlobalDiagonals();
    globalNumCols = input_matrix_->NumGlobalCols();
#endif // EEP
  }
  else //graph
  {
    myNZ = this->input_graph_->getLocalNumEntries();
    mySelfEdges = this->input_graph_->getLocalNumRows(); // Diagonals(); // EEP___
    globalNZ = this->input_graph_->getGlobalNumEntries();
    globalSelfEdges = this->input_graph_->getGlobalNumRows(); // Diagonals(); EEP___
    globalNumCols = this->input_graph_->getGlobalNumCols();
  }
  if (this->costs_.get() != 0)
  {
    throw std::runtime_error("EEP Nonzero costs_ pointer in isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::computeCost()");
#if 0 // EEP___
    numMyVWeights = this->costs_->getNumVertices();

    if (this->costs_->haveGraphEdgeWeights()){
        for (int i=0; i<numMyVWeights; i++){
          int gid = this->input_map_->getGlobalElement(i); // EEP
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
#endif // EEP___
  }
  else{
    this->costs_ = Teuchos::rcp(new CostDescriber<LocalOrdinal, GlobalOrdinal, Node>());
  }
  //comm.SumAll(&err, &gerr ,1); // EEP___

  if (gerr > 0){
    throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
  }

  int lval[4], gval[4];
  lval[0] = numMyVWeights;
  lval[1] = numMyGWeights;
  lval[2] = numMyHGWeights;

  //comm.SumAll(lval, gval, 3); // EEP___
  gval[0] = lval[0];
  gval[1] = lval[1];
  gval[2] = lval[2];

  int numVWeights = gval[0];
  int numGWeights = gval[1];
  int numHGWeights = gval[2];

  if ((numVWeights > 0) && (numVWeights != globalNumRows)){
    str2 = "Number of vertex weights supplied by application != number of rows";
    throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
  }
  if ((numGWeights > 0) && (numGWeights != (globalNZ - globalSelfEdges))){
    str2 = "Number of graph edge weights supplied by application != number of edges";
    throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
  }
  if ((numHGWeights > 0) && (numHGWeights < globalNumCols)){
    str2 = "Number of hyperedge weights supplied by application < number of columns";
    throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
  }

  this->costs_->setNumGlobalVertexWeights(numVWeights);
  this->costs_->setNumGlobalGraphEdgeWeights(numGWeights);
  this->costs_->setNumGlobalHypergraphEdgeWeights(numHGWeights);
  std::cout << "EEP Leaving isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::computeCost()" << std::endl;
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
void ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::preCheckPartition()
{
  std::string str1("Isorropia::ZoltanLibClass::precheckPartition ");
  std::string str2;

  //const Epetra_Comm &comm = input_map_->Comm();
  int localProc = this->input_map_->getComm()->getRank();
  int nprocs = this->input_map_->getComm()->getSize();

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
  int numrows = this->input_map_->getLocalNumElements(); // Elements(); // EEP___

  int myParts[] = {myGparts, myLparts};
  int maxParts[2];

  //comm.MaxAll(myParts, maxParts, 2); // EEP___
  maxParts[0] = myParts[0];
  maxParts[1] = myParts[1];

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

  if (maxLparts > 0) {
    //comm.SumAll(&myLparts, &sumLparts, 1); // EEP___
    sumLparts = myLparts;
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
        throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
      }

      if ((sumLparts > 0) && (sumLparts != maxGparts)){
        // This is an error because local partitions must sum to number of global
        str2 = "NUM_GLOBAL_PARTS not equal to sum of NUM_LOCAL_PARTS";
        throw std::runtime_error/*Isorropia::Exception*/(str1+str2);
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

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::
repartition(Teuchos::ParameterList& zoltanParamList,
            std::vector<int>& properties,
            int& exportsSize,
            std::vector<int>& imports)
{
  std::cout << "EEP Entering isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition()" << std::endl;
  zoltanParamList_  = zoltanParamList;

  // Avoid to construct import list.
  // Perhaps "PARTITION ASSIGNMENTS" will be better in term of performance.
  zoltanParamList_.set("RETURN_LISTS", "EXPORT AND IMPORT");

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 001" << std::endl;
  preCheckPartition(); // EEP check
  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 002" << std::endl;
  precompute();
  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 003"
            << ": this->numPartSizes = " << this->numPartSizes
            << std::endl;

  // Set part sizes

  if (this->numPartSizes > 0){
    int err;
    int *wgtIdx = NULL;

    err = zz_->LB_Set_Part_Sizes(1, this->numPartSizes, this->partGIDs, wgtIdx, this->partSizes);

    if (err != ZOLTAN_OK){
      throw std::runtime_error/*Isorropia::Exception*/("Error in LB_Set_Part_Sizes");
      return -1;
    }
  }

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 004" << std::endl;

  //Generate Load Balance
  int changes=0, num_gid_entries=0, num_lid_entries=0, num_import=0, num_export=0;
  ZOLTAN_ID_PTR import_global_ids=NULL, import_local_ids=NULL;
  ZOLTAN_ID_PTR export_global_ids=NULL, export_local_ids=NULL;
  int * import_procs=NULL, * export_procs=NULL;
  int *import_to_part=NULL, *export_to_part=NULL;

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 005" << std::endl;

  int err = zz_->LB_Partition( changes
                             , num_gid_entries
                             , num_lid_entries
                             , num_import
                             , import_global_ids
                             , import_local_ids
                             , import_procs
                             , import_to_part
                             , num_export
                             , export_global_ids
                             , export_local_ids
                             , export_procs
                             , export_to_part
                             );

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 006"
            << ": err = " << err
	    << std::endl;

  if (err != ZOLTAN_OK){
    throw std::runtime_error/*Isorropia::Exception*/("Error computing partitioning with Zoltan");
    return -1;
  }

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 007" << std::endl;

  exportsSize = num_export;
  imports.clear();
  imports.assign(import_global_ids, import_global_ids + num_import);

  properties.assign(num_obj_, queryObject_->RowMap().getComm()->getRank());

  for( int i = 0; i < num_export; ++i )
  {
    properties[export_local_ids[i]] = export_to_part[i];
  }

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 008" << std::endl;
  
  //Free Zoltan Data
  zz_->LB_Free_Part(&import_global_ids, &import_local_ids,
                     &import_procs, &import_to_part);
  zz_->LB_Free_Part(&export_global_ids, &export_local_ids,
                     &export_procs, &export_to_part);

  std::cout << "EEP In isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition(), pos 009" << std::endl;

  postcompute();

  std::cout << "EEP Leaving isorropia/src/tpetra/Isorropia_TpetraZoltanLib.hpp ZoltanLibClass<>::repartition()" << std::endl;
  return (0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::
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
      throw std::runtime_error/*Isorropia::Exception*/("Out of memory in ZoltanLibClass::color()");
      return -1;
    }

    int *intIds = nullptr; // queryObject_->RowMap().getLocalElementList(); // MyGlobalElements(); // EEP___
    for (int i=0; i < num_obj_; i++){
      gids[i] = (ZOLTAN_ID_TYPE)intIds[i];
    }
  }
  else{
    gids = nullptr; // (ZOLTAN_ID_PTR)queryObject_->RowMap().getLocalElementList(); // MyGlobalElements(); // EEP___
  }

  int err = zz_->Color(num_gid_entries, num_obj_, gids, &properties[0]);

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    delete [] gids;
  }

  if (err != ZOLTAN_OK){
    throw std::runtime_error/*Isorropia::Exception*/("Error computing coloring with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::
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
      throw std::runtime_error/*Isorropia::Exception*/("Out of memory in ZoltanLibClass::order()");
      return -1;
    }

    int *intIds = nullptr; // queryObject_->RowMap().MyGlobalElements(); // EEP___
    for (int i=0; i < num_obj_; i++){
      gids[i] = (ZOLTAN_ID_TYPE)intIds[i];
    }
  }
  else{
    gids = nullptr; // (ZOLTAN_ID_PTR)queryObject_->RowMap().MyGlobalElements(); // EEP___
  }

  int err = zz_->Order(num_gid_entries, num_obj_, gids, &properties[0], NULL);

  if (sizeof(ZOLTAN_ID_TYPE) != sizeof(int)){
    delete [] gids;
  }

  if (err != ZOLTAN_OK){
    throw std::runtime_error/*Isorropia::Exception*/("Error computing ordering with Zoltan");
    return -1;
  }

  postcompute();
  return (0);
}

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
int ZoltanLibClass<LocalOrdinal, GlobalOrdinal, Node>::postcompute()
{
  std::cout << "EEP Entering ZoltanLibClass<>::postcompute()" << std::endl;
  if (zz_)
    delete zz_;
  zz_ = NULL;

  std::cout << "EEP Leaving ZoltanLibClass<>::postcompute()" << std::endl;
  return (0);
}

}//namespace Tpetra
}//namespace Isorropia

#endif


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

