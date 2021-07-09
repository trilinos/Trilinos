// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file AdapterForTests.hpp
 *  \brief Generate Adapter for testing purposes.
 */


#ifndef ADAPTERFORTESTS
#define ADAPTERFORTESTS

#include <Zoltan2_Parameters.hpp>
#include <UserInputForTests.hpp>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_EvaluatePartition.hpp>

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_BasicVectorAdapter.hpp>

#ifdef HAVE_ZOLTAN2_PAMGEN
#include <Zoltan2_PamgenMeshAdapter.hpp>
#endif

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <string>
#include <iostream>
#include <vector>

using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::Comm;
using Teuchos::rcp;
using Teuchos::arcp;
using Teuchos::rcp_const_cast;
using Teuchos::ParameterList;
using std::string;
using namespace Zoltan2_TestingFramework;

// helper struct to store both an adapter and the coordinate adapter
struct AdapterWithTemplateName
{
  Zoltan2::BaseAdapterRoot * adapter = nullptr; // generic base class
  EAdapterType adapterType; // convert back to proper adapter type
};

struct AdapterWithOptionalCoordinateAdapter
{
  AdapterWithTemplateName main; // the main adapter - never null
  AdapterWithTemplateName coordinate; // can be null
};

/* \brief A class for constructing Zoltan2 input adapters */
class AdapterFactory{
public:
  /*! \brief A class method for constructing an input adapter
   *   defind in a parameter list.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the paramter list defining the data type to be used and
   *                 adapter type to be costructed
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterFactory(
    UserInputForTests *uinput, const ParameterList &pList,
    const RCP<const Comm<int> > &comm);

  ~AdapterFactory(); // handles deleting BaseAdapterRoot * data for adapter

  Zoltan2::BaseAdapterRoot * getMainAdapter() const {
    return adaptersSet.main.adapter;
  }

  EAdapterType getMainAdapterType() const {
    return adaptersSet.main.adapterType;
  }

  Zoltan2::BaseAdapterRoot * getCoordinateAdapter() const {
    return adaptersSet.coordinate.adapter;
  }

  EAdapterType getCoordinateAdapterType() const {
    return adaptersSet.coordinate.adapterType;
  }

private:
  AdapterWithOptionalCoordinateAdapter adaptersSet;

  /*! \brief Method to choose and call the correct constructor
   *   for a BasicIdentifierAdapter from a UserInputForTests input file.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithTemplateName
  getBasicIdentiferAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /*! \brief Method to choose and call the correct constructor
   *   for an Xpetra multi-vector adapter.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithTemplateName
  getXpetraMVAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /*! \brief Method to choose and call the correct constructor
   *   for an Xpetra crs graph adapter.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithOptionalCoordinateAdapter
  getXpetraCrsGraphAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /*! \brief Method to choose and call the correct constructor
   *   for an Xpetra Crs matrix adapter.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithOptionalCoordinateAdapter
  getXpetraCrsMatrixAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /*! \brief Method to choose and call the correct constructor
   *   for a basic vector adapter.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithTemplateName
  getBasicVectorAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);
  
  /*! \brief Method to choose and call the correct constructor
   *   for a Pamgen mesh adapter.
   *   \param[in] uinput is the data source for adapter
   *   \param[in] pList is the teuchos input parameter list
   *   \param[in] comm is the process communicator
   *
   * \return Ptr to the constructed adapter cast to the base class
   */
  AdapterWithTemplateName
  getPamgenMeshAdapterForInput(UserInputForTests *uinput,
    const ParameterList &pList, const RCP<const Comm<int> > &comm);

  /*! \brief Method to set up strided vector data from a multi-vector
   *  \param[in] data is the multi-vector
   *  \param[out] coords is the vector of strided coordinate data
   *  \param[out] strides is the vector of strides
   *  \param[in] stride is the stride to apply to data set
   *
   * \return
   */
  template <typename T>
  void InitializeVectorData(const RCP<T> &data,
                                   std::vector<const zscalar_t *> &coords,
                                   std::vector<int> & strides,
                                   int stride);
  
#ifdef HAVE_EPETRA_DATA_TYPES
  /*! \brief Method to set up strided vector data from a multi-vector
   *  \param[in] data is the epetra multi-vector
   *  \param[out] coords is the vector of strided coordinate data
   *  \param[out] strides is the vector of strides
   *  \param[in] stride is the stride to apply to data set
   *
   * \return
   */
  template <typename T>
  void InitializeEpetraVectorData(const RCP<T> &data,
                                         std::vector<const zscalar_t *> &coords,
                                         std::vector<int> & strides,
                                         int stride);
#endif
};


AdapterFactory::AdapterFactory(
  UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  if(!pList.isParameter("input adapter"))
  {
    std::cerr << "Input adapter unspecified" << std::endl;
    return;
  }

  // pick method for chosen adapter
  std::string input_adapter_name = pList.get<string>("input adapter");

  if(input_adapter_name == "BasicIdentifier")
    adaptersSet.main = getBasicIdentiferAdapterForInput(uinput, pList, comm);
  else if(input_adapter_name == "XpetraMultiVector")
    adaptersSet.main = getXpetraMVAdapterForInput(uinput, pList, comm);
  else if(input_adapter_name == "XpetraCrsGraph")
    adaptersSet = getXpetraCrsGraphAdapterForInput(uinput,pList, comm);
  else if(input_adapter_name == "XpetraCrsMatrix")
    adaptersSet = getXpetraCrsMatrixAdapterForInput(uinput,pList, comm);
  else if(input_adapter_name == "BasicVector")
    adaptersSet.main = getBasicVectorAdapterForInput(uinput,pList, comm);
  else if(input_adapter_name == "PamgenMesh")
    adaptersSet.main = getPamgenMeshAdapterForInput(uinput,pList, comm);

  if(adaptersSet.main.adapter == nullptr) {
    throw std::logic_error("AdapterFactory failed to create adapter!");
  }
}

AdapterFactory::~AdapterFactory() {
  if( adaptersSet.main.adapter ) {
    delete adaptersSet.main.adapter;
  }

  if( adaptersSet.coordinate.adapter ) {
    delete adaptersSet.coordinate.adapter;
  }
}


AdapterWithTemplateName
  AdapterFactory::getBasicIdentiferAdapterForInput(UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  AdapterWithTemplateName result;

  if(!pList.isParameter("data type"))
  {
    std::cerr << "Input data type unspecified" << std::endl;
    return result;
  }
  
  string input_type = pList.get<string>("data type"); // get the input type
  
  if (!uinput->hasInputDataType(input_type))
  {
    std::cerr << "Input type: " + input_type + " unavailable or misspelled."
              << std::endl; // bad type
    return result;
  }
  
  std::vector<const zscalar_t *> weights;
  std::vector<int> weightStrides;
  const zgno_t *globalIds = NULL;
  size_t localCount = 0;
  
  // get weights if any
  // get weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    // copy to weight
    size_t cols = vtx_weights->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
      weightStrides.push_back((int)vtx_weights->getStride());
    }
  }
  
  if(input_type == "coordinates")
  {
    RCP<tMVector_t> data = uinput->getUICoordinates();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "tpetra_vector")
  {
    RCP<tVector_t> data = uinput->getUITpetraVector();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "tpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<tMVector_t> data = uinput->getUITpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "tpetra_crs_graph")
  {
    RCP<tcrsGraph_t> data = uinput->getUITpetraCrsGraph();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "tpetra_crs_matrix")
  {
    RCP<tcrsMatrix_t> data = uinput->getUITpetraCrsMatrix();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "xpetra_vector")
  {
    RCP<xVector_t> data = uinput->getUIXpetraVector();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "xpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "xpetra_crs_graph")
  {
    RCP<xcrsGraph_t> data = uinput->getUIXpetraCrsGraph();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "xpetra_crs_matrix")
  {
    RCP<xcrsMatrix_t> data = uinput->getUIXpetraCrsMatrix();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_vector")
  {
    RCP<Epetra_Vector> data = uinput->getUIEpetraVector();
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = data->MyLength();
  }
  else if(input_type == "epetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = data->MyLength();
  }
  else if(input_type == "epetra_crs_graph")
  {
    RCP<Epetra_CrsGraph> data = uinput->getUIEpetraCrsGraph();
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = data->NumMyCols();
  }
  else if(input_type == "epetra_crs_matrix")
  {
    RCP<Epetra_CrsMatrix> data = uinput->getUIEpetraCrsMatrix();
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = data->NumMyCols();
  }
#endif

  result.adapterType = AT_basic_id_t;
  result.adapter = new Zoltan2_TestingFramework::basic_id_t(zlno_t(localCount),
                                                   globalIds,
                                                   weights,weightStrides);
  return result;
}


AdapterWithTemplateName AdapterFactory::getXpetraMVAdapterForInput(
  UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  AdapterWithTemplateName result;

  if(!pList.isParameter("data type"))
  {
    std::cerr << "Input data type unspecified" << std::endl;
    return result;
  }
  
  string input_type = pList.get<string>("data type");
  if (!uinput->hasInputDataType(input_type))
  {
    std::cerr << "Input type:" + input_type + ", unavailable or misspelled."
              << std::endl; // bad type
    return result;
  }
  
  std::vector<const zscalar_t *> weights;
  std::vector<int> weightStrides;
  
  // get weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    // copy to weight
    size_t weightsPerRow = vtx_weights->getNumVectors();
    for (size_t i = 0; i< weightsPerRow; i++) {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
      weightStrides.push_back(1);
    }
  }
  
  // set adapter
  if(input_type == "coordinates")
  {
    RCP<tMVector_t> data = uinput->getUICoordinates();
    RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
    if(weights.empty())
      result.adapter = new xMV_tMV_t(const_data);
    else {
      result.adapter = new xMV_tMV_t(const_data,weights,weightStrides);
    }
    result.adapterType = AT_xMV_tMV_t;
  }
  else if(input_type == "tpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<tMVector_t> data = uinput->getUITpetraMultiVector(nvec);
    RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
    if(weights.empty())
      result.adapter = new xMV_tMV_t(const_data);
    else
      result.adapter = new xMV_tMV_t(const_data,weights,weightStrides);
    result.adapterType = AT_xMV_tMV_t;
  }
  else if(input_type == "xpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(nvec);
    RCP<const xMVector_t> const_data = rcp_const_cast<const xMVector_t>(data);
    if(weights.empty())
      result.adapter = new xMV_xMV_t(const_data);
    else{
      result.adapter = new xMV_xMV_t(const_data,weights,weightStrides);
    }
    result.adapterType = AT_xMV_xMV_t;
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(nvec);
    RCP<const Epetra_MultiVector> const_data = rcp_const_cast<const Epetra_MultiVector>(data);
    
    if(weights.empty())
      result.adapter = new xMV_eMV_t(const_data);
    else
      result.adapter = new xMV_eMV_t(const_data,weights,weightStrides);
    result.adapterType = AT_xMV_eMV_t;
  }
#endif
  
  if(result.adapter == nullptr)
    std::cerr << "Input data chosen not compatible with xpetra multi-vector adapter." << std::endl;

  return result;
}


AdapterWithOptionalCoordinateAdapter AdapterFactory::getXpetraCrsGraphAdapterForInput(
  UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  
  AdapterWithOptionalCoordinateAdapter adapters;

  if(!pList.isParameter("data type"))
  {
    std::cerr << "Input data type unspecified" << std::endl;
    return adapters;
  }

  string input_type = pList.get<string>("data type");
  if (!uinput->hasInputDataType(input_type))
  {
    std::cerr << "Input type: " + input_type + ", unavailable or misspelled." 
              << std::endl; // bad type
    return adapters;
  }
  
  std::vector<const zscalar_t *> vtx_weights;
  std::vector<const zscalar_t *> edge_weights;
  std::vector<int> vtx_weightStride;
  std::vector<int> edge_weightStride;

  // get vtx weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights_tmp = uinput->getUIWeights();
    // copy to weight
    size_t weightsPerRow = vtx_weights_tmp->getNumVectors();
    for (size_t i = 0; i< weightsPerRow; i++) {
      vtx_weights.push_back(vtx_weights_tmp->getData(i).getRawPtr());
      vtx_weightStride.push_back(1);
    }
  }

  // get edge weights if any
  if(uinput->hasUIEdgeWeights())
  {
    RCP<tMVector_t> edge_weights_tmp = uinput->getUIEdgeWeights();
    // copy to weight
    size_t weightsPerRow = edge_weights_tmp->getNumVectors();
    for (size_t i = 0; i< weightsPerRow; i++) {
      edge_weights.push_back(edge_weights_tmp->getData(i).getRawPtr());
      edge_weightStride.push_back(1);
    }
  }
  
  // make the coordinate adapter
  // get an adapter for the coordinates
  // need to make a copy of the plist and change the vector type
  Teuchos::ParameterList pCopy(pList);
  pCopy = pCopy.set<std::string>("data type","coordinates");

  // for coordinate adapter
  #define SET_COORDS_INPUT_1(adapterClass)                                     \
        auto * ca = dynamic_cast<adapterClass*>(adapters.coordinate.adapter);  \
        if(!ca) {throw std::logic_error( "Coordinate adapter case failed!" );} \
        ia->setCoordinateInput(ca);

  if(input_type == "tpetra_crs_graph")
  {
    RCP<tcrsGraph_t> data = uinput->getUITpetraCrsGraph();
    RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(data);

    xCG_tCG_t * ia = new xCG_tCG_t(const_data,(int)vtx_weights.size(),(int)edge_weights.size());
    adapters.main.adapterType = AT_xCG_tCG_t;
    adapters.main.adapter = ia;

    if(!vtx_weights.empty()) {
      for(int i = 0; i < (int)vtx_weights.size(); i++)
        ia->setVertexWeights(vtx_weights[i],vtx_weightStride[i],i);
    }
    
    if(!edge_weights.empty()) {
      for(int i = 0; i < (int)edge_weights.size(); i++)
        ia->setEdgeWeights(edge_weights[i],edge_weightStride[i],i);
    }
    
    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_1);
    }
  }
  else if(input_type == "xpetra_crs_graph")
  {
    RCP<xcrsGraph_t> data = uinput->getUIXpetraCrsGraph();
    RCP<const xcrsGraph_t> const_data = rcp_const_cast<const xcrsGraph_t>(data);

    xCG_xCG_t * ia = new xCG_xCG_t(const_data, (int)vtx_weights.size(), (int)edge_weights.size());
    adapters.main.adapterType = AT_xCG_xCG_t;
    adapters.main.adapter = ia;
    if(!vtx_weights.empty())
    {
      for(int i = 0; i < (int)vtx_weights.size(); i++)
        ia->setVertexWeights(vtx_weights[i],vtx_weightStride[i],i);
    }
    
    if(!edge_weights.empty())
    {
      for(int i = 0; i < (int)edge_weights.size(); i++)
        ia->setEdgeWeights(edge_weights[i],edge_weightStride[i],i);
    }

    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_1);
    }
  }
#ifdef HAVE_EPETRA_DATA_TYPES

  else if(input_type == "epetra_crs_graph")
  {
    RCP<Epetra_CrsGraph> data = uinput->getUIEpetraCrsGraph();
    RCP<const Epetra_CrsGraph> const_data = rcp_const_cast<const Epetra_CrsGraph>(data);
    xCG_eCG_t * ia = new xCG_eCG_t(const_data,(int)vtx_weights.size(),(int)edge_weights.size());
    adapters.main.adapterType = AT_xCG_eCG_t;
    adapters.main.adapter = ia;
    if(!vtx_weights.empty())
    {
      for(int i = 0; i < (int)vtx_weights.size(); i++)
        ia->setVertexWeights(vtx_weights[i],vtx_weightStride[i],i);
    }
    
    if(!edge_weights.empty())
    {
      for(int i = 0; i < (int)edge_weights.size(); i++)
        ia->setEdgeWeights(edge_weights[i],edge_weightStride[i],i);
    }

    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_1);
    }
  }
#endif
  
  if(adapters.main.adapter == nullptr) {
    std::cerr << "Input data chosen not compatible with "
              << "XpetraCrsGraph adapter." << std::endl;
    return adapters;
  }

  return adapters;
}


AdapterWithOptionalCoordinateAdapter AdapterFactory::getXpetraCrsMatrixAdapterForInput(
  UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  AdapterWithOptionalCoordinateAdapter adapters;

  if(!pList.isParameter("data type"))
  {
    std::cerr << "Input data type unspecified" << std::endl;
    return adapters;
  }
  
  string input_type = pList.get<string>("data type");
  if (!uinput->hasInputDataType(input_type))
  {
    std::cerr << "Input type:" + input_type + ", unavailable or misspelled."
              << std::endl; // bad type
    return adapters;
  }
  
  std::vector<const zscalar_t *> weights;
  std::vector<int> strides;
  
  // get weights if any
  if(uinput->hasUIWeights())
  {
    if(comm->getRank() == 0) std::cout << "Have weights...." << std::endl;
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    
    // copy to weight
    int weightsPerRow = (int)vtx_weights->getNumVectors();
    for (int i = 0; i< weightsPerRow; i++)
    {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
      strides.push_back(1);
    }
    
  }

  // make the coordinate adapter
  // get an adapter for the coordinates
  // need to make a copy of the plist and change the vector type
  Teuchos::ParameterList pCopy(pList);
  pCopy = pCopy.set<std::string>("data type","coordinates");

  // for coordinate adapter
  #define SET_COORDS_INPUT_2(adapterClass)                                     \
        auto * ca = dynamic_cast<adapterClass*>(adapters.coordinate.adapter);  \
        if(!ca) {throw std::logic_error( "Coordinate adapter case failed!" );} \
        ia->setCoordinateInput(ca);

  // set adapter
  if(input_type == "tpetra_crs_matrix")
  {
    if(comm->getRank() == 0) std::cout << "Make tpetra crs matrix adapter...." << std::endl;
    
    // get pointer to data
    RCP<tcrsMatrix_t> data = uinput->getUITpetraCrsMatrix();
    RCP<const tcrsMatrix_t> const_data = rcp_const_cast<const tcrsMatrix_t>(data); // const cast data
    
    // new adapter
    xCM_tCM_t *ia = new xCM_tCM_t(const_data, (int)weights.size());
    adapters.main.adapterType = AT_xCM_tCM_t;
    adapters.main.adapter = ia;

    // if we have weights set them
    if(!weights.empty())
    {
      for(int i = 0; i < (int)weights.size(); i++)
        ia->setWeights(weights[i],strides[i],i);
    }

    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_2);
    }
  }
  else if(input_type == "xpetra_crs_matrix")
  {
    RCP<xcrsMatrix_t> data = uinput->getUIXpetraCrsMatrix();
    RCP<const xcrsMatrix_t> const_data = rcp_const_cast<const xcrsMatrix_t>(data);
    
    // new adapter
    xCM_xCM_t *ia = new xCM_xCM_t(const_data, (int)weights.size());
    adapters.main.adapterType = AT_xCM_xCM_t;
    adapters.main.adapter = ia;

    // if we have weights set them
    if(!weights.empty())
    {
      for(int i = 0; i < (int)weights.size(); i++)
         ia->setWeights(weights[i],strides[i],i);
    }

    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_2);
    }
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_crs_matrix")
  {
    RCP<Epetra_CrsMatrix> data = uinput->getUIEpetraCrsMatrix();
    RCP<const Epetra_CrsMatrix> const_data = rcp_const_cast<const Epetra_CrsMatrix>(data);
    
    // new adapter
    xCM_eCM_t *ia = new xCM_eCM_t(const_data, (int)weights.size());
    adapters.main.adapterType = AT_xCM_eCM_t;
    adapters.main.adapter = ia;

    // if we have weights set them
    if(!weights.empty())
    {
      for(int i = 0; i < (int)weights.size(); i++)
         ia->setWeights(weights[i],strides[i],i);
    }

    if (uinput->hasUICoordinates()) {
      adapters.coordinate = getXpetraMVAdapterForInput(uinput, pCopy, comm);
      Z2_TEST_UPCAST_COORDS(adapters.coordinate.adapterType, SET_COORDS_INPUT_2);
    }
  }
#endif
  
  if(adapters.main.adapter == nullptr)
  {
    std::cerr << "Input data chosen not compatible with "
              << "XpetraCrsMatrix adapter." << std::endl;
    return adapters;
  }

  return adapters;
}

AdapterWithTemplateName AdapterFactory::getBasicVectorAdapterForInput(
  UserInputForTests *uinput,
  const ParameterList &pList,
  const RCP<const Comm<int> > &comm)
{
  
  AdapterWithTemplateName result;

  if(!pList.isParameter("data type"))
  {
    std::cerr << "Input data type unspecified" << std::endl;
    return result;
  }
  
  string input_type = pList.get<string>("data type");
  if (!uinput->hasInputDataType(input_type))
  {
    std::cerr << "Input type:" + input_type + ", unavailable or misspelled."
              << std::endl; // bad type
    return result;
  }
  
  std::vector<const zscalar_t *> weights;
  std::vector<int> weightStrides;
  const zgno_t * globalIds;
  zlno_t localCount = 0;
  
  // get weights if any
  // get weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    // copy to weight
    size_t cols = vtx_weights->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
      weightStrides.push_back(1);
    }
  }
  
  // get vector stride
  int stride = 1;
  if(pList.isParameter("stride"))
    stride = pList.get<int>("stride");

  result.adapterType = AT_basic_vector_adapter;

  if(input_type == "coordinates")
  {
    RCP<tMVector_t> data = uinput->getUICoordinates();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = static_cast<zlno_t>(data->getLocalLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeVectorData(data,coords,entry_strides,stride);
    

    
    if (weights.empty()) {
      size_t dim = coords.size(); //BDD add NULL for constructor call
      size_t push_null = 3-dim;
      for (size_t i = 0; i < push_null; i ++) coords.push_back(NULL);
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(
                                                     zlno_t(localCount),
                                                     globalIds,
                                                     coords[0],
                                                     coords[1],coords[2],
                                                     stride, stride, stride);
    } else if (weights.size() == 1) {
      size_t dim = coords.size(); //BDD add NULL for constructor call
      size_t push_null = 3-dim;
      for (size_t i = 0; i < push_null; i ++) coords.push_back(NULL);
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(
                                                     zlno_t(localCount),
                                                     globalIds,
                                                     coords[0],
                                                     coords[1],coords[2],
                                                     stride, stride, stride,
                                                     true,
                                                     weights[0],
                                                     weightStrides[0]);
    } else { // More than one weight per ID
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(
                                                     zlno_t(localCount),
                                                     globalIds,
                                                     coords, entry_strides,
                                                     weights, weightStrides);
    }
  }
  else if(input_type == "tpetra_vector")
  {
    RCP<tVector_t> data = uinput->getUITpetraVector();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = static_cast<zlno_t>(data->getLocalLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeVectorData(data,coords,entry_strides,stride);
    
    if(weights.empty())
    {
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0]);
    }else{
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0],
                                                     true,
                                                     weights[0],
                                                     weightStrides[0]);
      
    }
    
  }
  else if(input_type == "tpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    
    RCP<tMVector_t> data = uinput->getUITpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = static_cast<zlno_t>(data->getLocalLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeVectorData(data,coords,entry_strides,stride);
    
    result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                   coords, entry_strides,
                                                   weights,weightStrides);
    
  }
  else if(input_type == "xpetra_vector")
  {
    RCP<xVector_t> data = uinput->getUIXpetraVector();
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = static_cast<zlno_t>(data->getLocalLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeVectorData(data,coords,entry_strides,stride);
    
    if(weights.empty())
    {
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0]);
    }else{
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0],
                                                     true,
                                                     weights[0],
                                                     weightStrides[0]);
      
    }
  }
  else if(input_type == "xpetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = static_cast<zlno_t>(data->getLocalLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeVectorData(data,coords,entry_strides,stride);
    if(comm->getRank() == 0) std::cout << "size of entry strides: " << entry_strides.size() << std::endl;
    if(comm->getRank() == 0) std::cout << "size of coords: " << coords.size() << std::endl;
    
    // make vector!
    result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                   coords, entry_strides,
                                                   weights,weightStrides);
  }
  
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_vector")
  {
    RCP<Epetra_Vector> data = uinput->getUIEpetraVector();
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = static_cast<zlno_t>(data->MyLength());
    
    // get strided data
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeEpetraVectorData(data,coords,entry_strides,stride);
    if(weights.empty())
    {
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0]);
    }else{
      result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                     coords[0], entry_strides[0],
                                                     true,
                                                     weights[0],
                                                     weightStrides[0]);
      
    }
    
    //    delete [] epetravectors;
  }
  else if(input_type == "epetra_multivector")
  {
    int nvec = pList.get<int>("vector_dimension");
    RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(nvec);
    globalIds = (zgno_t *)data->Map().MyGlobalElements();
    localCount = data->MyLength();
    
    std::vector<const zscalar_t *> coords;
    std::vector<int> entry_strides;
    InitializeEpetraVectorData(data,coords,entry_strides,stride);
    
    // make vector!
    result.adapter = new Zoltan2_TestingFramework::basic_vector_adapter(localCount, globalIds,
                                                   coords, entry_strides,
                                                   weights,weightStrides);
  }
  
#endif
  
  return result;
}

template <typename T>
void AdapterFactory::InitializeVectorData(const RCP<T> &data,
                                           std::vector<const zscalar_t *> &coords,
                                           std::vector<int> & strides,
                                           int stride)
{
  // set up adapter data
  size_t localCount = data->getLocalLength();
  size_t nvecs = data->getNumVectors();
  size_t vecsize = data->getNumVectors() * data->getLocalLength();
//    printf("Number of vectors by data: %zu\n", nvecs);
  //  printf("Size of data: %zu\n", vecsize);
  
  ArrayRCP<zscalar_t> *petravectors = new ArrayRCP<zscalar_t>[nvecs];
  
  //  printf("Getting t-petra vectors...\n");
  for (size_t i = 0; i < nvecs; i++)
    petravectors[i] = data->getDataNonConst(i);
  
  // debugging
  //  for (size_t i = 0; i < nvecs; i++){
  //    printf("Tpetra vector %zu: {",i);
  //
  //    for (size_t j = 0; j < localCount; j++)
  //    {
  //      printf("%1.2g ",petravectors[i][j]);
  //    }
  //    printf("}\n");
  //  }
  
  size_t idx = 0;
  zscalar_t *coordarr = new zscalar_t[vecsize];
  
  if(stride == 1 || stride != (int)nvecs)
  {
    for (size_t i = 0; i < nvecs; i++) {
      for (size_t j = 0; j < localCount; j++) {
        coordarr[idx++] = petravectors[i][j];
      }
    }
  }else
  {
    for (size_t j = 0; j < localCount; j++) {
      for (size_t i = 0; i < nvecs; i++) {
        coordarr[idx++] = petravectors[i][j];
      }
    }
  }
  
  // debugging
  //  printf("Made coordarr : {");
  //  for (zlno_t i = 0; i < vecsize; i++){
  //    printf("%1.2g ",coordarr[i]);
  //  }
  //  printf("}\n");
  
  // always build for dim 3
  coords = std::vector<const zscalar_t *>(nvecs);
  strides = std::vector<int>(nvecs);
  
  for (size_t i = 0; i < nvecs; i++) {
    if(stride == 1)
      coords[i] = &coordarr[i*localCount];
    else
      coords[i] = &coordarr[i];
    
    strides[i] = stride;
  }
  
  // debugging
  //  printf("Made coords...\n");
  //  for (size_t i = 0; i < nvecs; i++){
  //    const zscalar_t * tmp = coords[i];
  //    printf("coord %zu: {",i);
  //    for(size_t j = 0; j < localCount; j++)
  //    {
  //      printf("%1.2g ", tmp[j]);
  //    }
  //    printf("}\n");
  //  }
  
  //  printf("clean up coordarr and tpetravectors...\n\n\n");
  delete [] petravectors;
}

#ifdef HAVE_EPETRA_DATA_TYPES

template <typename T>
void AdapterFactory::InitializeEpetraVectorData(const RCP<T> &data,
                                                 std::vector<const zscalar_t *> &coords,
                                                 std::vector<int> & strides,
                                                 int stride){
  size_t localCount = data->MyLength();
  size_t nvecs = data->NumVectors();
  size_t vecsize = nvecs * localCount;
  
  //  printf("Number of vectors by data: %zu\n", nvecs);
  //  printf("Size of data: %zu\n", vecsize);
  
  std::vector<zscalar_t *> epetravectors(nvecs);
  zscalar_t ** arr;
  //  printf("get data from epetra vector..\n");
  data->ExtractView(&arr);
  
  for(size_t k = 0; k < nvecs; k++)
  {
    epetravectors[k] = arr[k];
  }
  
  size_t idx = 0;
  basic_vector_adapter::scalar_t *coordarr =
    new basic_vector_adapter::scalar_t[vecsize];

  if(stride == 1 || stride != (int)nvecs)
  {
    for (size_t i = 0; i < nvecs; i++) {
      for (size_t j = 0; j < localCount; j++) {
        coordarr[idx++] = epetravectors[i][j];
      }
    }
  }else
  {
    for (size_t j = 0; j < localCount; j++) {
      for (size_t i = 0; i < nvecs; i++) {
        coordarr[idx++] = epetravectors[i][j];
      }
    }
  }
  
  // debugging
//  printf("Made coordarr : {");
//  for (zlno_t i = 0; i < vecsize; i++){
//    printf("%1.2g ",coordarr[i]);
//  }
//  printf("}\n");
  
  coords = std::vector<const zscalar_t *>(nvecs);
  strides = std::vector<int>(nvecs);
  
  for (size_t i = 0; i < nvecs; i++) {
    if(stride == 1)
      coords[i] = &coordarr[i*localCount];
    else
      coords[i] = &coordarr[i];
    
    strides[i] = stride;
  }
  
//  printf("Made coords...\n");
//  for (size_t i = 0; i < nvecs; i++){
//    const zscalar_t * tmp = coords[i];
//    printf("coord %zu: {",i);
//    for(size_t j = 0; j < localCount; j++)
//    {
//      printf("%1.2g ", tmp[j]);
//    }
//    printf("}\n");
//  }
  
}
#endif


// pamgen adapter
AdapterWithTemplateName
AdapterFactory::getPamgenMeshAdapterForInput(UserInputForTests *uinput,
                                              const ParameterList &pList,
                                              const RCP<const Comm<int> > &comm)
{
  AdapterWithTemplateName result;

#ifdef HAVE_ZOLTAN2_PAMGEN
  if(uinput->hasPamgenMesh())
  {
    if(uinput->hasPamgenMesh())
    {
//      if(comm->getRank() == 0) std::cout << "Have pamgen mesh, constructing adapter...." << std::endl;
      result.adapter =
        new pamgen_adapter_t(*(comm.get()), "region");
      result.adapterType = AT_pamgen_adapter_t;
//      if(comm->getRank() == 0)
//        ia->print(0);
    }
  }else{
    std::cerr << "Pamgen mesh is unavailable for PamgenMeshAdapter!"
              << std::endl;
  }
  
  return result;
#else
  throw std::runtime_error("Pamgen input requested but Trilinos is not "
                           "built with Pamgen");
#endif
}
#endif


