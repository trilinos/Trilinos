// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
//  Zoltan2_TestInterface.h
//  Zoltan2TestDriver
//
//  Created by Bradley Davidson on 7/6/15.
//  Copyright (c) 2015 TXCorp. All rights reserved.
//

#ifndef Zoltan2TestDriver_Zoltan2_TestInterface_h
#define Zoltan2TestDriver_Zoltan2_TestInterface_h

#include <Zoltan2_config.h>
#include <Zoltan2_TestHelpers.hpp>
#include <UserInputForTests.hpp>
#include <Teuchos_DefaultComm.hpp>

/*! \brief Zoltan2_TestInteface defines the api for all Tests.
 
 A test derived from Zoltan2Test can be easily
 called by the Zoltan2 test driver.
 
 */
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ParameterList;

class Zoltan2Test{
    
public:
    
    /*! \brief Destructor
     */
    virtual ~Zoltan2Test() {};
    
    /*! \breif Run the test
     \param \c uinput user input helper object
     \param \c communicator object
     */
    virtual void Run(const ParameterList &params,const RCP<const Teuchos::Comm<int> > & comm) = 0;
    
    /*! \brief Did pass?
     */
    virtual bool didPass() = 0;
    
};

#endif
