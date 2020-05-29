//
//  Zoltan2_Tests.h
//  Zoltan2TestDriver
//
//  Created by Bradley Davidson on 7/6/15.
//  Copyright (c) 2015 TXCorp. All rights reserved.
//

#ifndef Zoltan2TestDriver_Zoltan2_Tests_h
#define Zoltan2TestDriver_Zoltan2_Tests_h

#include "Zoltan2_MeshCoordinateTest.hpp"
#include <string>
using std::string;
Zoltan2Test * getZoltan2Test(const string &testType)
{
    if(testType.find("MeshCoordinateTest") != string::npos)
    {
        return new MeshCoordinateTest();
    }
    
    return nullptr;
}
#endif
