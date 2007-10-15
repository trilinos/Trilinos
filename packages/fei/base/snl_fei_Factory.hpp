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
#include <snl_fei_Solver.hpp>
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
    fei::SharedPtr<fei::Factory> clone() const
      {
	fei::SharedPtr<fei::Factory> factory;
	if (wrapper_.get() != NULL) {
	  factory.reset(new snl_fei::Factory(comm_, wrapper_));
	}
	else if (lsc_.get() != NULL) {
	  factory.reset(new snl_fei::Factory(comm_, lsc_));
	}
	else if (feData_.get() != NULL) {
	  factory.reset(new snl_fei::Factory(comm_, feData_, nodeIDType_));
	}

	return(factory);
      }

    /** Implementation of fei::Factory::parameters() */
    virtual void parameters(const fei::ParameterSet& parameterset)
      {
        fei::Factory::parameters(parameterset);

	int err = 0;
	if (lsc_.get() != NULL || feData_.get() != NULL) {
	  int numParams = 0;
	  const char** paramStrings = NULL;
	  std::vector<std::string> stdstrings;
	  fei::utils::convert_ParameterSet_to_strings(&parameterset, stdstrings);
	  fei::utils::strings_to_char_ptrs(stdstrings, numParams, paramStrings);

	  char** nc_paramStrings = const_cast<char**>(paramStrings);
	  if (lsc_.get() != NULL) {
	    err += lsc_->parameters(numParams, nc_paramStrings);
	  }
	  if (feData_.get() != NULL) {
	    err += feData_->parameters(numParams, nc_paramStrings);
	  }

	  delete [] paramStrings;

	  if (err != 0) {
	    FEI_OSTRINGSTREAM osstr;
	    osstr << "snl_fei::Factory::parameters received err="<<err
		  << " from either feiData_->parameters or lsc_->parameters.";
	    throw fei::Exception(osstr.str());
	  }
	}

	parameterset.getIntParamValue("outputLevel", outputLevel_);
      }

    /** Implementation of fei::MatrixGraph::Factory::createMatrixGraph.
    */
    virtual fei::SharedPtr<fei::MatrixGraph>
      createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                        fei::SharedPtr<fei::VectorSpace> columnSpace,
                        const char* name)
      {
        static fei::MatrixGraph_Impl2::Factory factory;
        matrixGraph_ = factory.createMatrixGraph(rowSpace,
                                              columnSpace, name);
        return(matrixGraph_);
      }

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector>
      createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
		   int numVectors=1)
      {
	(void)vecSpace;
	fei::SharedPtr<fei::Vector> dummy;
	if (matrixGraph_.get() == NULL) {
	  FEI_CERR << "snl_fei::Factory ERROR: when using LinearSystemCore or FiniteElementData"
	       << ", you must create a MatrixGraph before you can create vectors"<<FEI_ENDL;
	  return(dummy);
	}

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	return( broker_->createVector() );
      }

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector> createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
						     bool isSolutionVector,
						     int numVectors=1)
      {
	fei::SharedPtr<fei::Vector> dummy;
	(void)vecSpace;
	if (matrixGraph_.get() == NULL) {
	  FEI_CERR << "snl_fei::Factory ERROR: when using LinearSystemCore"
	       << ", you must create a MatrixGraph before you can create vectors"<<FEI_ENDL;
	  return(dummy);
	}

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	return( broker_->createVector(isSolutionVector) );
      }

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector> createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
						     int numVectors=1)
      {
	matrixGraph_ = matrixGraph;

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	return( broker_->createVector() );
      }

    /** Implementation of fei::Vector::Factory::createVector() */
    virtual fei::SharedPtr<fei::Vector> createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
						     bool isSolutionVector,
						     int numVectors=1)
      {
	matrixGraph_ = matrixGraph;

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	return( broker_->createVector(isSolutionVector) );
      }

    /** Implementation of fei::Matrix::Factory::createMatrix() */
    virtual fei::SharedPtr<fei::Matrix> createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
      {
	matrixGraph_ = matrixGraph;
	fei::SharedPtr<fei::Matrix> mptr;

	if (matrixGraph_.get() == NULL) {
	  FEI_CERR << "snl_fei::Factory ERROR: when using LinearSystemCore"
	       << ", you must create a MatrixGraph before you can create matrices"<<FEI_ENDL;
	  return(mptr);
	}

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	broker_->setMatrixGraph(matrixGraph);

	return(broker_->createMatrix());
      }

    /** Implementation of 
	fei::LinearSystem::Factory::createLinearSystem() */
    virtual fei::SharedPtr<fei::LinearSystem>
      createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
      {
	matrixGraph_ = matrixGraph;

	if (matrixGraph_.get() == NULL) {
	  FEI_CERR << "snl_fei::Factory ERROR: you may not create a LinearSystem with "
	       << "a NULL MatrixGraph object." << FEI_ENDL;
	  fei::SharedPtr<fei::LinearSystem> linsys;
	  return(linsys);
	}

        if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
            reducer_.get() == NULL) {
          reducer_ = matrixGraph_->getReducer();
        }

	createBroker(matrixGraph_);

	broker_->setMatrixGraph(matrixGraph);

	return( broker_->createLinearSystem() );
      }

    /** Implementation of fei::Solver::Factory::createSolver() */
    virtual fei::SharedPtr<fei::Solver> createSolver(const char* name=0)
      {
	fei::SharedPtr<fei::Solver> solver(new snl_fei::Solver);
	return(solver);
      }

    /** get LibraryWrapper attribute (power-users only) */
    fei::SharedPtr<LibraryWrapper> get_LibraryWrapper() const
      {
	return( wrapper_ );
      }

    int getOutputLevel()
      {
	return(outputLevel_);
      }

  private:
    int createBroker(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
    {
      int err = -1;
      if (lsc_.get() != NULL) {
	err = createBroker_LinSysCore(matrixGraph, lsc_);
      }
      if (feData_.get() != NULL) {
	err = createBroker_FEData(matrixGraph, feData_);
      }

      return(err);
    }

    int createBroker_LinSysCore(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
				fei::SharedPtr<LinearSystemCore> lsc)
      {
	if (broker_.get() == NULL) {
	  fei::SharedPtr<snl_fei::Broker>
	  brokerptr(new snl_fei::Broker_LinSysCore(lsc, matrixGraph, reducer_));
	  broker_ = brokerptr;
	}

	return(0);
      }

    int createBroker_FEData(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
			    fei::SharedPtr<FiniteElementData> feData)
      {
	if (broker_.get() == NULL) {
	  fei::SharedPtr<snl_fei::Broker>
	    brokerptr(new snl_fei::Broker_FEData(feData, matrixGraph,
						 nodeIDType_));
	  broker_ = brokerptr;
	}

	return(0);
      }

    MPI_Comm comm_;
    fei::SharedPtr<snl_fei::Broker> broker_;
    fei::SharedPtr<fei::MatrixGraph> matrixGraph_;
    fei::SharedPtr<fei::Reducer> reducer_;

    int nodeIDType_;

    fei::SharedPtr<LinearSystemCore> lsc_;
    fei::SharedPtr<FiniteElementData> feData_;
    fei::SharedPtr<LibraryWrapper> wrapper_;
    int outputLevel_;
  };//class Factory
}//namespace snl_fei

#endif // _snl_fei_Factory_hpp_
