//
//  AdapterForTests.h
//  Zoltan2TestDriver
//
//  Created by Bradley Davidson on 7/10/15.
//  Copyright (c) 2015 TXCorp. All rights reserved.
//

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
#include <Zoltan2_BasicVectorAdapter.hpp>
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
    typedef Zoltan2::BasicIdentifierAdapter<userTypes_t> basic_adapter;
    typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> xpetra_mv_adapter;
    /*! \brief Method to choose and call the coorect constructor
     *      for a BasicIdentifierAdapter from a UserInputForTests input file.
     *   \param uinput is the UserInputForTestsInputForTestObject
     *   \param  pList is the teuchos input parameter list
     *   \param  adapter is a reference to the input adapter to be constructed.
     */
    static basic_adapter getBasicIdentiferAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
    
    /*! \brief Method to choose and call the coorect constructor
     *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
     *   \param uinput is the UserInputForTestsInputForTestObject
     *   \param  pList is the teuchos input parameter list
     *   \param  adapter is a reference to the input adapter to be constructed.
     */
    static xpetra_mv_adapter getXpetraMVAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
    
    /*! \brief Method to choose and call the coorect constructor
     *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
     *   \param uinput is the UserInputForTestsInputForTestObject
     *   \param  pList is the teuchos input parameter list
     *   \param  adapter is a reference to the input adapter to be constructed.
     */
    static void
    getXpetraCRSGraphAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
    
    /*! \brief Method to choose and call the coorect constructor
     *      for a XpetraMultiVectorAdapter from a UserInputForTests input file.
     *   \param uinput is the UserInputForTestsInputForTestObject
     *   \param  pList is the teuchos input parameter list
     *   \param  adapter is a reference to the input adapter to be constructed.
     */
    static void
    getXpetraCRSMatrixAdapterForInput(UserInputForTests *uinput, const ParameterList &pList);
};

AdapterForTests::basic_adapter AdapterForTests::getBasicIdentiferAdapterForInput(UserInputForTests *uinput, const ParameterList &pList)
{
    if(!pList.isParameter("inputType"))
        throw std::runtime_error("Input data type not specified");
    
    string input_type = pList.get<string>("inputType"); // get the input type
    
    if (!uinput->hasInputDataType(input_type))
        throw std::runtime_error("Input type not avaible, or misspelled."); // bad type
    
    vector<const zscalar_t *> weights;
    std::vector<int> weightStrides;
    const zzgid_t * globalIds;
    size_t localCount;
    
    // get weights if any
    if(uinput->hasUIWeights())
    {
        cout << "\t\tProblem has vtx weights! adding weights..." << endl;
        RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
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
        RCP<tMVector_t> data = uinput->getUITpetraMultiVector(1);
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
        RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(0);
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
        globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
        localCount = data->getLocalLength();
    }
    else if(input_type == "epetra_multivector")
    {
        RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(0);
        globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
        localCount = data->getLocalLength();
    }
    else if(input_type == "epetra_crs_graph")
    {
        RCP<Epetra_CrsGraph> data = uinput->getUIEpetraCrsGraph();
        globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
        localCount = data->getNodeNumCols();
    }
    else if(input_type == "epetra_crs_matrix")
    {
        RCP<Epetra_CrsMatrix> data = uinput->getUIEpetraCrsMatrix();
        globalIds = (zzgid_t *)data->getMap()->getNodeElementList().getRawPtr();
        localCount = data->getNodeNumCols();
    }
#endif
    
    return Zoltan2::BasicIdentifierAdapter<userTypes_t>(zlno_t(localCount),globalIds,weights,weightStrides);
}


AdapterForTests::xpetra_mv_adapter AdapterForTests::getXpetraMVAdapterForInput(UserInputForTests *uinput, const ParameterList &pList)
{
    
    if(!pList.isParameter("inputType"))
        throw std::runtime_error("Input data type not specified");
    
    string input_type = pList.get<string>("inputType");
    if (!uinput->hasInputDataType(input_type))
        throw std::runtime_error("Input type not avaible, or misspelled.");
    
    vector<const zscalar_t *> weights;
    std::vector<int> weightStrides;
    
    // get weights if any
    if(uinput->hasUIWeights())
    {
        cout << "\t\tProblem has vtx weights! adding weights..." << endl;
        RCP<tMVector_t> vtx_weights = uinput->getUIWeights();
    }
    
    // set adapter
        if(input_type == "coordinates")
        {
            RCP<tMVector_t> data = uinput->getUICoordinates();
            RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
            if(weights.empty())
                return Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data);
            else
                return Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data,weights,weightStrides);
        }
        else if(input_type == "tpetra_multivector")
        {
            RCP<tMVector_t> data = uinput->getUITpetraMultiVector(1);
            RCP<const tMVector_t> const_data = rcp_const_cast<const tMVector_t>(data);
            if(weights.empty())
                return Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data);
            else
                return Zoltan2::XpetraMultiVectorAdapter<tMVector_t>(const_data,weights,weightStrides);
        }
//        else if(input_type == "xpetra_multivector")
//        {
//            RCP<xMVector_t> data = uinput->getUIXpetraMultiVector(1);
//            cout << "\t\t...assigning pointer to XpetraMultiVectorAdapter...." << endl;
//            RCP<const xMVector_t> const_data = rcp_const_cast<const xMVector_t>(data);
//            if(weights.empty())
//                return Zoltan2::XpetraMultiVectorAdapter<xMVector_t>(const_data);
//            else
//                return Zoltan2::XpetraMultiVectorAdapter<xMVector_t>(const_data,weights,weightStrides);
//        }
#ifdef HAVE_EPETRA_DATA_TYPES
    
    else if(input_type == "epetra_multivector")
    {
        RCP<Epetra_MultiVector> data = uinput->getUIEpetraMultiVector(1);
        RCP<const Epetra_MultiVector> const_data = rcp_const_cast<const Epetra_MultiVector>(data);
        if(weights.empty())
            return Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector>(const_data);
        else
            return Zoltan2::XpetraMultiVectorAdapter<Epetra_MultiVector>(const_data,weights,weightStrides);
    }
#endif
    
    throw std::runtime_error("Input data chosen not compatible with xpetra multi-vector adapter.");
}

//
//void AdapterForTests::getXpetraCRSGraphAdapterForInput(UserInputForTests *uinput)
//{
//    cout << "CREATING XPETRA CRS GRAPH ADAPTER ...." <<endl;
//}
//
//
//void AdapterForTests::getXpetraCRSMatrixAdapterForInput(UserInputForTests *uinput,
//                                                                           const ParameterList &pList)
//{
//    cout << "CREATING XPETRA CRS MATRIX ADAPTER ...." <<endl;
//}

#endif
