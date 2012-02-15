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


#ifndef _snl_fei_Broker_hpp_
#define _snl_fei_Broker_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>

namespace fei {
  class VectorSpace;
  class MatrixGraph;
  class Vector;
  class Matrix;
  class LinearSystem;
}//namespace fei

namespace snl_fei {

  /** Internal interface. Similar to a factory, for creating
      Matrix/Vector/etc instances which are views of LinearSystemCore
      implementations. One of these will be instantiated and
      used internally by snl_fei::Factory.
  */
  class Broker {
  public:
    /** Usual virtual destructor. */
    virtual ~Broker(){}

    /** Produce an instance of an fei::Vector. This overloading of the
	create() method is for use by Broker implementations that are
	dispensing 'views' of vectors that reside in LinearSystemCore or
	FiniteElementData container implementations. In those cases, there is
	a distinction that must be made between solution-vectors and
	rhs-vectors.

	@param isSolutionVector
     */
    virtual fei::SharedPtr<fei::Vector> createVector(bool isSolutionVector=false) = 0;

    /** Produce an instance of an fei::Matrix
     */
    virtual fei::SharedPtr<fei::Matrix> createMatrix() = 0;

    /** Produce an instance of an fei::LinearSystem
     */
    virtual fei::SharedPtr<fei::LinearSystem> createLinearSystem() = 0;

    /** Set the MatrixGraph object used by this broker. */
    virtual void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph) = 0;
  };//class Broker
}//namespace snl_fei

#endif // _snl_fei_Broker_hpp_
