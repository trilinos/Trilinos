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
#include <Zoltan2_PartitioningSolutionQuality.hpp>

#include <Zoltan2_BasicIdentifierAdapter.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <string>
#include <iostream>
#include <vector>

using namespace std;
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


class AdapterForTests{
public:
  
  typedef UserInputForTests::tcrsMatrix_t tcrsMatrix_t;
  typedef UserInputForTests::tcrsGraph_t tcrsGraph_t;
  typedef UserInputForTests::tVector_t tVector_t;
  typedef UserInputForTests::tMVector_t tMVector_t;
  
  typedef UserInputForTests::xcrsMatrix_t xcrsMatrix_t;
  typedef UserInputForTests::xcrsGraph_t xcrsGraph_t;
  typedef UserInputForTests::xVector_t xVector_t;
  typedef UserInputForTests::xMVector_t xMVector_t;
  
  typedef Zoltan2::BasicUserTypes<zscalar_t, zzgid_t, zlno_t, zgno_t> userTypes_t;
  typedef Zoltan2::BaseAdapter<userTypes_t> base_adapter_t;
  typedef Zoltan2::BasicIdentifierAdapter<userTypes_t> basic_id_t;
  typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> xpetra_mv_adapter;
  
  static base_adapter_t* getAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
  
private:
  /*! \brief Method to choose and call the coorect constructor
   *      for a BasicIdentifierAdapter from a UserInputForTests input file.
   *   \param uinput is the UserInputForTestsInputForTestObject
   *   \param  pList is the teuchos input parameter list
   *   \param  adapter is a reference to the input adapter to be constructed.
   */
  static base_adapter_t*
  getBasicIdentiferAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
  
  /*! \brief Method to choose and call the coorect constructor
   *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
   *   \param uinput is the UserInputForTestsInputForTestObject
   *   \param  pList is the teuchos input parameter list
   *   \param  adapter is a reference to the input adapter to be constructed.
   */
  static base_adapter_t*
  getXpetraMVAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
  
  /*! \brief Method to choose and call the coorect constructor
   *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
   *   \param uinput is the UserInputForTestsInputForTestObject
   *   \param  pList is the teuchos input parameter list
   *   \param  adapter is a reference to the input adapter to be constructed.
   */
  static base_adapter_t*
  getXpetraCrsGraphAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
  
  /*! \brief Method to choose and call the coorect constructor
   *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
   *   \param uinput is the UserInputForTestsInputForTestObject
   *   \param  pList is the teuchos input parameter list
   *   \param  adapter is a reference to the input adapter to be constructed.
   */
  static base_adapter_t*
  getXpetraCrsMatrixAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
};


AdapterForTests::base_adapter_t * AdapterForTests::getAdapterForInput(UserInputForTests *uinput, const ParameterList &pList)
{
  // ask chaco gen what data types are available! -> for debugging
  
  if(uinput->hasUICoordinates())
    std::cout << "coordinates available" << std::endl;
  if(uinput->hasUIWeights())
    std::cout << "vtx weights available" << std::endl;
  if(uinput->hasUIEdgeWeights())
    std::cout << "edge weights available" << std::endl;
  if(uinput->hasUITpetraCrsMatrix())
    std::cout << "tpet crs matrix available" << std::endl;
  if(uinput->hasUITpetraCrsGraph())
    std::cout << "tpet crs graph available" << std::endl;
  if(uinput->hasUITpetraVector())
    std::cout << "teptra vector available" << std::endl;
  if(uinput->hasUITpetraMultiVector())
    std::cout << "tpet mv available" << std::endl;
  if(uinput->hasUIXpetraCrsMatrix())
    std::cout << "xpet crs matrix available" << std::endl;
  if(uinput->hasUIXpetraCrsGraph())
    std::cout << "xpet crs graph available" << std::endl;
  if(uinput->hasUIXpetraVector())
    std::cout << "xpet vector available" << std::endl;
  if(uinput->hasUIXpetraMultiVector())
    std::cout << "xpet mv available" << std::endl;
#ifdef HAVE_EPETRA_DATA_TYPES
  if(uinput->hasUIEpetraCrsGraph())
    std::cout << "epet crs graph available" << std::endl;
  if(uinput->hasUIEpetraCrsMatrix())
    std::cout << "epet crs matrix available" << std::endl;
  if(uinput->hasUIEpetraVector())
    std::cout << "epet vector available" << std::endl;
  if(uinput->hasUIEpetraMultiVector())
    std::cout << "epet mv available" << std::endl;
#endif
  
  
  if(!pList.isParameter("inputAdapter"))
    throw std::runtime_error("Input adapter not specified");
  
  // pick method for chosen adapter
  string adapter_name = pList.get<string>("inputAdapter");
  AdapterForTests::base_adapter_t * ia = nullptr; // input adapter
  if(adapter_name == "BasicIdentifier")
    ia = AdapterForTests::getBasicIdentiferAdapterForInput(uinput, pList);
  else if(adapter_name == "XpetraMultiVector")
    ia = AdapterForTests::getXpetraMVAdapterForInput(uinput, pList);
  //    else if(adapter_name == "XpetraCrsGraph")
  //        problemWithXpetraCrsGraphAdapter<inputAdapter_t, data_t>(ia,zoltan2params,comm);
  //    else if(adapter_name == "XpetraCrsMatrix")
  //        problemWithXpetraCrsMatrixAdapter<inputAdapter_t, data_t>(ia,zoltan2params,comm);
  else
    throw std::runtime_error("Input adapter type not avaible, or misspelled.");
  
  return ia;
}


AdapterForTests::base_adapter_t * AdapterForTests::getBasicIdentiferAdapterForInput(UserInputForTests *uinput, const ParameterList &pList)
{
  
  if(!pList.isParameter("inputType"))
    throw std::runtime_error("Input data type not specified");
  
  string input_type = pList.get<string>("inputType"); // get the input type
  
  if (!uinput->hasInputDataType(input_type))
    throw std::runtime_error("Input type not avaible, or misspelled."); // bad type
  
  vector<const zscalar_t *> weights;
  std::vector<int> weightStrides;
  const zzgid_t * globalIds;
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
    }
  }
  
  if(input_type == "coordinates")
  {
    RCP<tMVector_t> data = uinput->getUICoordinates();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
    cout <<"\t\tgot coordinates..." << endl;
  }
  else if(input_type == "tpetra_vector")
  {
    RCP<tVector_t> data = uinput->getUITpetraVector();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "tpetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<tMVector_t> data = uinput->getUITpetraMultiVector(nvec);
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "tpetra_crs_graph")
  {
    RCP<tcrsGraph_t> data = uinput->getUITpetraCrsGraph();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "tpetra_crs_matrix")
  {
    RCP<tcrsMatrix_t> data = uinput->getUITpetraCrsMatrix();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "xpetra_vector")
  {
    RCP<xVector_t> data = uinput->getUIXpetraVector();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "xpetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(nvec);
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getLocalLength();
  }
  else if(input_type == "xpetra_crs_graph")
  {
    RCP<xcrsGraph_t> data = uinput->getUIXpetraCrsGraph();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
  else if(input_type == "xpetra_crs_matrix")
  {
    RCP<xcrsMatrix_t> data = uinput->getUIXpetraCrsMatrix();
    globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
    localCount = data->getNodeNumCols();
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  else if(input_type == "epetra_vector")
  {
    RCP<Epetra_Vector> data = uinput->getUIEpetraVector();
    globalIds = (zzgid_t *)data->Map().MyGlobalElements();
    localCount = data->MyLength();
  }
  else if(input_type == "epetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(nvec);
    globalIds = (zzgid_t *)data->Map().MyGlobalElements();
    localCount = data->MyLength();
  }
  else if(input_type == "epetra_crs_graph")
  {
    RCP<Epetra_CrsGraph> data = uinput->getUIEpetraCrsGraph();
    globalIds = (zzgid_t *)data->Map().MyGlobalElements();
    localCount = data->NumMyCols();
  }
  else if(input_type == "epetra_crs_matrix")
  {
    RCP<Epetra_CrsMatrix> data = uinput->getUIEpetraCrsMatrix();
    globalIds = (zzgid_t *)data->Map().MyGlobalElements();
    localCount = data->NumMyCols();
  }
#endif
  
  if(localCount == 0) return nullptr;
  return reinterpret_cast<AdapterForTests::base_adapter_t *>( new AdapterForTests::basic_id_t(zlno_t(localCount),globalIds,weights,weightStrides));
}


AdapterForTests::base_adapter_t * AdapterForTests::getXpetraMVAdapterForInput(UserInputForTests *uinput, const ParameterList &pList)
{
  
  if(!pList.isParameter("inputType"))
    throw std::runtime_error("Input data type not specified");
  
  string input_type = pList.get<string>("inputType");
  if (!uinput->hasInputDataType(input_type))
    throw std::runtime_error("Input type not avaible, or misspelled.");
  
  AdapterForTests::base_adapter_t * adapter = nullptr;
  vector<const zscalar_t *> weights;
  std::vector<int> weightStrides;
  
  // get weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    
    // copy to weight
    size_t cols = vtx_weights->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
    }
  }
  
  // set adapter
  if(input_type == "coordinates")
  {
    RCP<tMVector_t> data = uinput->getUICoordinates();
    RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
    if(weights.empty())
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data));
    else
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data,weights,weightStrides));
  }
  else if(input_type == "tpetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<tMVector_t> data = uinput->getUITpetraMultiVector(nvec);
    RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
    if(weights.empty())
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data));
    else
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data,weights,weightStrides));
  }
  else if(input_type == "xpetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(nvec);
    RCP<const xMVector_t> const_data = rcp_const_cast<const xMVector_t>(data);
    if(weights.empty())
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<xMVector_t>(const_data));
    else{
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(new Zoltan2::XpetraMultiVectorAdapter<xMVector_t>(const_data,weights,weightStrides));
    }
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  
  else if(input_type == "epetra_multivector")
  {
    int nvec = pList.get<int>("number_of_vectors");
    RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(nvec);
    RCP<const Epetra_MultiVector> const_data = rcp_const_cast<const Epetra_MultiVector>(data);
    
    if(weights.empty())
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(
                                                                    new Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector>(const_data));
    else
      adapter = reinterpret_cast<AdapterForTests::base_adapter_t *>(
                                                                    new Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector>(const_data,weights,weightStrides));
  }
#endif
  
  if(adapter == nullptr)
    throw std::runtime_error("Input data chosen not compatible with xpetra multi-vector adapter.");
  else
    return adapter;
}


AdapterForTests::base_adapter_t * AdapterForTests::getXpetraCrsGraphAdapterForInput(UserInputForTests *uinput,
                                                                                    const ParameterList &pList)
{
  if(!pList.isParameter("inputType"))
    throw std::runtime_error("Input data type not specified");
  
  string input_type = pList.get<string>("inputType");
  if (!uinput->hasInputDataType(input_type))
    throw std::runtime_error("Input type not avaible, or misspelled.");
  
  
  AdapterForTests::base_adapter_t * adapter = nullptr;
  vector<const zscalar_t *> vtx_weights;
  vector<const zscalar_t *> edge_weights;
  int vtx_weightStride = 1;
  int edge_weightStride = 1;
  int vtx_idx = 0;
  int edge_idx = 0;
  
  // get vtx weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights_tmp = uinput->getUIWeights();
    
    // copy to weight
    size_t cols = vtx_weights_tmp->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      vtx_weights.push_back(vtx_weights_tmp->getData(i).getRawPtr());
    }
  }
  
  // get edge weights if any
  if(uinput->hasUIEdgeWeights())
  {
    RCP<tMVector_t> edge_weights_tmp = uinput->getUIEdgeWeights();
    
    // copy to weight
    size_t cols = edge_weights_tmp->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      edge_weights.push_back(edge_weights_tmp->getData(i).getRawPtr());
    }
  }
  
  
  // set adapter
  if(input_type == "tpetra_crs_graph")
  {
    typedef Zoltan2::XpetraCrsGraphAdapter<tcrsGraph_t, tcrsGraph_t> problem_t;
    
    RCP<tcrsGraph_t> data = uinput->getUITpetraCrsGraph();
    RCP<const tcrsGraph_t> const_data = rcp_const_cast<const tcrsGraph_t>(data);
    problem_t *ia = new problem_t(const_data,(int)vtx_weights.size(),(int)edge_weights.size());
    
    if(!vtx_weights.empty()) ia->setVertexWeights(vtx_weights[0],vtx_weightStride,vtx_idx);
    if(!edge_weights.empty()) ia->setEdgeWeights(edge_weights[0],edge_weightStride,edge_idx);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
  }
  else if(input_type == "xpetra_crs_graph")
  {
    typedef Zoltan2::XpetraCrsGraphAdapter<xcrsGraph_t, xcrsGraph_t> problem_t;
    
    RCP<xcrsGraph_t> data = uinput->getUIXpetraCrsGraph();
    RCP<const xcrsGraph_t> const_data = rcp_const_cast<const xcrsGraph_t>(data);
    problem_t *ia = new problem_t(const_data, (int)vtx_weights.size(), (int)edge_weights.size());
    
    if(!vtx_weights.empty()) ia->setVertexWeights(vtx_weights[0],vtx_weightStride,vtx_idx);
    if(!edge_weights.empty()) ia->setEdgeWeights(edge_weights[0],edge_weightStride,edge_idx);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  
  else if(input_type == "epetra_crs_graph")
  {
    typedef Zoltan2::XpetraCrsGraphAdapter<Epetra_CrsGraph, Epetra_CrsGraph> problem_t;
    
    RCP<Epetra_CrsGraph> data = uinput->getUIEpetraCrsGraph();
    RCP<const Epetra_CrsGraph> const_data = rcp_const_cast<const Epetra_CrsGraph>(data);
    problem_t *ia = new problem_t(const_data,(int)vtx_weights.size(),(int)edge_weights.size());
    
    if(!vtx_weights.empty()) ia->setVertexWeights(vtx_weights[0],vtx_weightStride,vtx_idx);
    if(!edge_weights.empty()) ia->setEdgeWeights(edge_weights[0],edge_weightStride,edge_idx);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
    
  }
#endif
  
  if(adapter == nullptr)
    throw std::runtime_error("Input data chosen not compatible with xpetra multi-vector adapter.");
  else
    return adapter;
  
}


AdapterForTests::base_adapter_t * AdapterForTests::getXpetraCrsMatrixAdapterForInput(UserInputForTests *uinput,
                                                                                     const ParameterList &pList)
{
  if(!pList.isParameter("inputType"))
    throw std::runtime_error("Input data type not specified");
  
  string input_type = pList.get<string>("inputType");
  if (!uinput->hasInputDataType(input_type))
    throw std::runtime_error("Input type not avaible, or misspelled.");
  
  
  AdapterForTests::base_adapter_t * adapter = nullptr;
  vector<const zscalar_t *> weights;
  int weightStride = 1;
  
  // get weights if any
  if(uinput->hasUIWeights())
  {
    RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    
    // copy to weight
    size_t cols = vtx_weights->getNumVectors();
    for (size_t i = 0; i< cols; i++) {
      weights.push_back(vtx_weights->getData(i).getRawPtr());
    }
  }
  
  // set adapter
  if(input_type == "tpetra_crs_matrix")
  {
    typedef Zoltan2::XpetraCrsMatrixAdapter<tcrsMatrix_t, tcrsMatrix_t> problem_t;
    
    RCP<tcrsMatrix_t> data = uinput->getUITpetraCrsMatrix();
    RCP<const tcrsMatrix_t> const_data = rcp_const_cast<const tcrsMatrix_t>(data);
    problem_t *ia = new problem_t(const_data, (int)weights.size());
    if(!weights.empty()) ia->setWeights(weights[0],weightStride);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
  }
  else if(input_type == "xpetra_crs_matrix")
  {
    typedef Zoltan2::XpetraCrsMatrixAdapter<xcrsMatrix_t, xcrsMatrix_t> problem_t;
    
    RCP<xcrsMatrix_t> data = uinput->getUIXpetraCrsMatrix();
    RCP<const xcrsMatrix_t> const_data = rcp_const_cast<const xcrsMatrix_t>(data);
    problem_t *ia = new problem_t(const_data, (int)weights.size());
    if(!weights.empty()) ia->setWeights(weights[0],weightStride);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
  }
#ifdef HAVE_EPETRA_DATA_TYPES
  
  else if(input_type == "epetra_crs_matrix")
  {
    typedef Zoltan2::XpetraCrsMatrixAdapter<Epetra_CrsMatrix, Epetra_CrsMatrix> problem_t;
    
    RCP<Epetra_CrsMatrix> data = uinput->getUIEpetraCrsMatrix();
    RCP<const Epetra_CrsMatrix> const_data = rcp_const_cast<const Epetra_CrsMatrix>(data);
    problem_t *ia = new problem_t(const_data, (int)weights.size());
    if(!weights.empty()) ia->setWeights(weights[0],weightStride);
    
    adapter =  reinterpret_cast<AdapterForTests::base_adapter_t *>(ia);
  }
#endif
  
  if(adapter == nullptr)
    throw std::runtime_error("Input data chosen not compatible with xpetra multi-vector adapter.");
  else
    return adapter;
}

#endif
