/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

