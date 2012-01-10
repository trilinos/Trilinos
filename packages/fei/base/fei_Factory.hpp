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


#ifndef _fei_Factory_hpp_
#define _fei_Factory_hpp_

#include "fei_macros.hpp"
#include "fei_mpi.h"
#include "fei_VectorSpace.hpp"
#include "fei_MatrixGraph.hpp"
#include "fei_Matrix.hpp"
#include "fei_Vector.hpp"
#include "fei_LinearSystem.hpp"
#include "fei_Solver.hpp"
#include "fei_LibraryWrapper.hpp"
#include "FEI.hpp"

namespace fei {
  //first, a forward declaration...
  class ParameterSet;

  /** Interface for creating fei:: instances.
      In all cases, input arguments (arguments required to construct the
      requested class) are followed by the result or output argument.

      This interface inherits the various fei:: factory interfaces as a
      convenience mechanism, so that user code can deal with one factory
      object instead of a different factory for each class. In addition
      to inheriting the fei:: factories, this interface also provides
      methods for creating instances of the 'old' FEI class.
   */
  class Factory : public virtual fei::VectorSpace::Factory,
                  public virtual fei::MatrixGraph::Factory,
                  public virtual fei::Matrix::Factory,
                  public virtual fei::Vector::Factory,
                  public virtual fei::LinearSystem::Factory,
                  public virtual fei::Solver::Factory {
  public:
    /** constructor */
    Factory(MPI_Comm comm);

    /** virtual destructor */
    virtual ~Factory();

    /** Create and return a new Factory of the same type. */
    virtual fei::SharedPtr<Factory> clone() const = 0;

    /** Set parameters.
     */
    virtual void parameters(const fei::ParameterSet& paramset);

    /** Produce an instance of the "old" FEI class (implements the FEI 2.1
	interface specification).

        This function is virtual, but not pure-virtual. An implementation
        is provided by this class, and can be inherited by derived classes
        if desired.
    */
    virtual fei::SharedPtr<FEI> createFEI(fei::SharedPtr<LibraryWrapper> wrapper,
					  MPI_Comm comm);

    /** Produce an instance of the "old" FEI class (implements the FEI 2.1
	interface specification).

        This function is virtual, but not pure-virtual. An implementation
        is provided by this class, and can be inherited by derived classes
        if desired.
    */
    virtual fei::SharedPtr<FEI> createFEI(MPI_Comm comm);

    /** Query screen output-level (set by parameter-string "outputLevel n"
	via parameters())
     */
    virtual int getOutputLevel() const = 0;

   private:
    Factory();
    Factory(const Factory& src);
  };//class Factory
}//namespace fei

#endif // _fei_Factory_hpp_
