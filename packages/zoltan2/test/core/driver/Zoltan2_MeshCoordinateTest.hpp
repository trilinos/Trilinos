// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
//  Zoltan2_MeshCoordinateTest.h
//  Zoltan2TestDriver
//
//  Created by Bradley Davidson on 7/6/15.
//  Copyright (c) 2015 TXCorp. All rights reserved.
//

#ifndef Zoltan2TestDriver_Zoltan2_MeshCoordinateTest_h
#define Zoltan2TestDriver_Zoltan2_MeshCoordinateTest_h

#include "Zoltan2_TestInterface.hpp"
#include <Zoltan2_BasicVectorAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_PartitioningProblem.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

typedef Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t> tMVector_t;
typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> myTypes_t;


class MeshCoordinateTest : public Zoltan2Test{
    
public:
    /*! \brief Default Constructor
     */
    MeshCoordinateTest();
    
    /*! \brief Destructor
     */
    ~MeshCoordinateTest() {};
    
    /*! \breif Run the test
     \param \c uinput user input helper object
     \param \c communicator object
     */
    void Run(const ParameterList &params,const RCP<const Teuchos::Comm<int> > & comm);
    
    /*! \brief Did pass?
     */
    bool didPass();
    
private:
    bool success;
};


MeshCoordinateTest::MeshCoordinateTest(){
    this->success = false;
};

void MeshCoordinateTest::Run(const ParameterList &params,
                             const RCP<const Teuchos::Comm<int> > & comm)
{
    const ParameterList &input = params.sublist("TestParameters");
    
    UserInputForTests uinput(input,comm,true, true);
    if(!uinput.hasUICoordinates()) return;
    
    RCP<tMVector_t> coords = uinput.getUICoordinates();
    
    size_t localCount = coords->getLocalLength();
    
    zscalar_t *x=NULL, *y=NULL, *z=NULL;
    x = coords->getDataNonConst(0).getRawPtr();
    y = coords->getDataNonConst(1).getRawPtr();
    z = coords->getDataNonConst(2).getRawPtr();
    
    const zgno_t *globalIds = coords->getMap()->getLocalElementList().getRawPtr();
    typedef Zoltan2::BasicVectorAdapter<tMVector_t> inputAdapter_t;
    
    inputAdapter_t ia(localCount, globalIds, x, y, z, 1, 1, 1);
    
//    ParameterList zoltan2params(params.sublist("Zoltan2Parameters"));
    const ParameterList &zoltan2params = params.sublist("Zoltan2Parameters");
#ifdef HAVE_ZOLTAN2_MPI
    Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, const_cast<ParameterList *>(&zoltan2params), MPI_COMM_WORLD);
#else
    Zoltan2::PartitioningProblem<inputAdapter_t> problem(&ia, const_cast<ParameterList *>(&zoltan2params));
#endif
    
    problem.solve();
    
    this->success = true;
    
}

bool MeshCoordinateTest::didPass(){return this->success;}


#endif
