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


#ifndef _snl_fei_Factory_hpp_
#define _snl_fei_Factory_hpp_

#include <fei_macros.hpp>
#include <fei_Factory.hpp>
#include <fei_VectorSpace.hpp>
#include <snl_fei_Broker_FEData.hpp>
#include <snl_fei_Broker_LinSysCore.hpp>
#include <fei_Solver.hpp>
#include <fei_LinearSystemCore.hpp>
#include <fei_LibraryWrapper.hpp>
#include <fei_utils.hpp>
#include <fei_Reducer.hpp>
#include <fei_MatrixReducer.hpp>
#include <fei_VectorReducer.hpp>
#include <fei_ParameterSet.hpp>
#include <fei_MatrixGraph_Impl2.hpp>

#undef fei_file
#define fei_file "snl_fei_Factory.hpp"
#include <fei_ErrMacros.hpp>

namespace snl_fei {

  /** snl_fei:: implementation of the various fei:: Factory interfaces.
   */
  class Factory : public virtual fei::Factory {
  public:
    /** Constructor */
    Factory(MPI_Comm comm, fei::SharedPtr<LibraryWrapper> wrapper);

    /** Constructor */
    Factory(MPI_Comm comm, fei::SharedPtr<LinearSystemCore> lsc);

    /** Constructor */
    Factory(MPI_Comm comm,
            fei::SharedPtr<FiniteElementData> feData, int nodeIDType);

    /** Destructor */
    virtual ~Factory();


    /** Implementation of fei::Factory::clone() */
    fei::SharedPtr<fei::Factory> clone() const;

    /** Implementation of fei::Factory::parameters() */
    virtual void parameters(const fei::ParameterSet& parameterset);

    /** Implementation of fei::MatrixGraph::Factory::createMatrixGraph.
    */
    virtual fei::SharedPtr<fei::MatrixGraph>
      createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                        fei::SharedPtr<fei::VectorSpace> columnSpace,
                        const char* name);

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector>
      createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
		   int numVectors=1);

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector> createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
						     bool isSolutionVector,
						     int numVectors=1);

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector>
      createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                   int numVectors=1);

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector>
      createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
		     bool isSolutionVector,
		     int numVectors=1);

    /** Implementation of fei::Matrix::Factory::createMatrix() */
    virtual fei::SharedPtr<fei::Matrix> createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

    /** Implementation of 
	fei::LinearSystem::Factory::createLinearSystem() */
    virtual fei::SharedPtr<fei::LinearSystem>
      createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph);

    /** Implementation of fei::Solver::Factory::createSolver() */
    virtual fei::SharedPtr<fei::Solver> createSolver(const char* name=0);

    /** get LibraryWrapper attribute (power-users only) */
    fei::SharedPtr<LibraryWrapper> get_LibraryWrapper() const;

    int getOutputLevel() const;

  private:
    int createBroker(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

    int createBroker_LinSysCore(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
				fei::SharedPtr<LinearSystemCore> lsc);

    int createBroker_FEData(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
			    fei::SharedPtr<FiniteElementData> feData);

    MPI_Comm comm_;
    fei::SharedPtr<snl_fei::Broker> broker_;
    fei::SharedPtr<fei::MatrixGraph> matrixGraph_;
    fei::SharedPtr<fei::Reducer> reducer_;

    int nodeIDType_;

    fei::SharedPtr<LinearSystemCore> lsc_;
    fei::SharedPtr<FiniteElementData> feData_;
    fei::SharedPtr<LibraryWrapper> wrapper_;
    int outputLevel_;
    bool blockMatrix_;
  };//class Factory
}//namespace snl_fei

#endif // _snl_fei_Factory_hpp_

