/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_constants_hpp_
#define _fei_constants_hpp_

/*
   In this file we define some constants to use as parameters to
   some fei functions.
   These constants are primarily used as 'fieldType' arguments to
   fei::VectorSpace::defineFields and FEI::initFields.
   If defining a vector-field for displacement, use DISPLACEMENT. If
   defining separate scalar fields for the components of displacement,
   then use DISPLACEMENT_X, etc.
   Most of the names below are self-explanatory. PRESSURE refers to either
   a vector-field for pressure, a nodal pressure variable, or the constant
   coefficient for a pressure field that is discontinuous in each element.
   PRESSURE_X, PRESSURE_Y, and PRESSURE_Z refer to the X, Y, and Z coefficients
   for a linearly varying pressure field defined separately in each element.
*/

namespace fei {

const int DISPLACEMENT     =  0;
const int DISPLACEMENT_X   =  0;
const int DISPLACEMENT_Y   =  1;
const int DISPLACEMENT_Z   =  2;
const int ROTATION         =  3;
const int ROTATION_X       =  3;
const int ROTATION_Y       =  4;
const int ROTATION_Z       =  5;
const int VELOCITY         =  6;
const int VELOCITY_X       =  6;
const int VELOCITY_Y       =  7;
const int VELOCITY_Z       =  8;
const int PRESSURE         =  9;
const int PRESSURE_X       = 10;
const int PRESSURE_Y       = 11;
const int PRESSURE_Z       = 12;
const int TEMPERATURE      = 13;

const int UNKNOWN          = 20;

}//namespace fei

#endif

