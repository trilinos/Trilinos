/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef _snl_fei_Factory_cpp_
#define _snl_fei_Factory_cpp_

#include <fei_macros.hpp>

#include <snl_fei_Factory.hpp>

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<LibraryWrapper> wrapper)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(0),
    lsc_(),
    feData_(),
    wrapper_(wrapper),
    outputLevel_(0),
    blockMatrix_(false)
{
  if (wrapper_.get() != NULL) {
    lsc_ = wrapper->getLinearSystemCore();
    feData_ = wrapper->getFiniteElementData();
  }
}

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<LinearSystemCore> lsc)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(0),
    lsc_(lsc),
    feData_(),
    wrapper_(),
    outputLevel_(0),
    blockMatrix_(false)
{
}

//----------------------------------------------------------------------------
snl_fei::Factory::Factory(MPI_Comm comm,
                          fei::SharedPtr<FiniteElementData> feData, int nodeIDType)
  : fei::Factory(comm),
    comm_(comm),
    broker_(),
    matrixGraph_(),
    nodeIDType_(nodeIDType),
    lsc_(),
    feData_(feData),
    wrapper_(NULL),
    outputLevel_(0),
    blockMatrix_(false)
{
}

//----------------------------------------------------------------------------
snl_fei::Factory::~Factory()
{
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Factory>
snl_fei::Factory::clone() const
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

//----------------------------------------------------------------------------
void
snl_fei::Factory::parameters(const fei::ParameterSet& parameterset)
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
      throw std::runtime_error(osstr.str());
    }
  }

  parameterset.getIntParamValue("outputLevel", outputLevel_);

  const fei::Param* param = 0;
  fei::Param::ParamType ptype = fei::Param::BAD_TYPE;

  param = parameterset.get("BLOCK_GRAPH");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    blockMatrix_ = true;
  }

  param = parameterset.get("BLOCK_MATRIX");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    if (ptype == fei::Param::BOOL) {
      blockMatrix_ = param->getBoolValue();
    }
    else {
      blockMatrix_ = true;
    }
  }

  param = parameterset.get("AZ_matrix_type");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype != fei::Param::BAD_TYPE) {
    if (ptype == fei::Param::STRING) {
      if (param->getStringValue() == "AZ_VBR_MATRIX") {
        blockMatrix_ = true;
      }
    }
  }
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::MatrixGraph>
snl_fei::Factory::createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                        fei::SharedPtr<fei::VectorSpace> columnSpace,
                        const char* name)
{
  static fei::MatrixGraph_Impl2::Factory factory;
  matrixGraph_ = factory.createMatrixGraph(rowSpace, columnSpace, name);
  return(matrixGraph_);
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Vector>
snl_fei::Factory::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                   int numVectors)
{
  (void)vecSpace;
  fei::SharedPtr<fei::Vector> dummy;
  if (matrixGraph_.get() == NULL) {
    fei::console_out() << "snl_fei::Factory ERROR: when using LinearSystemCore or FiniteElementData"
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

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Vector>
snl_fei::Factory::createVector(fei::SharedPtr<fei::VectorSpace> vecSpace,
                                                     bool isSolutionVector,
                                                     int numVectors)
{
  fei::SharedPtr<fei::Vector> dummy;
  (void)vecSpace;
  if (matrixGraph_.get() == NULL) {
    fei::console_out() << "snl_fei::Factory ERROR: when using LinearSystemCore"
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

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Vector>
snl_fei::Factory::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                                     int numVectors)
{
  matrixGraph_ = matrixGraph;

  if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
      reducer_.get() == NULL) {
    reducer_ = matrixGraph_->getReducer();
  }

  createBroker(matrixGraph_);

  return( broker_->createVector() );
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Vector>
snl_fei::Factory::createVector(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                     bool isSolutionVector,
                     int numVectors)
{
  matrixGraph_ = matrixGraph;

  if (matrixGraph_->getGlobalNumSlaveConstraints() > 0 &&
      reducer_.get() == NULL) {
    reducer_ = matrixGraph_->getReducer();
  }

  createBroker(matrixGraph_);

  return( broker_->createVector(isSolutionVector) );
}

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Matrix>
snl_fei::Factory::createMatrix(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  matrixGraph_ = matrixGraph;
  fei::SharedPtr<fei::Matrix> mptr;

  if (matrixGraph_.get() == NULL) {
    fei::console_out() << "snl_fei::Factory ERROR: when using LinearSystemCore"
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

//----------------------------------------------------------------------------
fei::SharedPtr<fei::LinearSystem>
snl_fei::Factory::createLinearSystem(fei::SharedPtr<fei::MatrixGraph>& matrixGraph)
{
  matrixGraph_ = matrixGraph;

  if (matrixGraph_.get() == NULL) {
    fei::console_out() << "snl_fei::Factory ERROR: you may not create a LinearSystem with "
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

//----------------------------------------------------------------------------
fei::SharedPtr<fei::Solver>
snl_fei::Factory::createSolver(const char* name)
{
  fei::SharedPtr<fei::Solver> solver(new fei::Solver);
  return(solver);
}

//----------------------------------------------------------------------------
fei::SharedPtr<LibraryWrapper>
snl_fei::Factory::get_LibraryWrapper() const
{
  return( wrapper_ );
}

//----------------------------------------------------------------------------
int snl_fei::Factory::getOutputLevel() const
{
  return(outputLevel_);
}

//----------------------------------------------------------------------------
int
snl_fei::Factory::createBroker(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
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

//----------------------------------------------------------------------------
int
snl_fei::Factory::createBroker_LinSysCore(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                                fei::SharedPtr<LinearSystemCore> lsc)
{
  if (broker_.get() == NULL) {
    fei::SharedPtr<snl_fei::Broker> brokerptr(new snl_fei::Broker_LinSysCore(lsc, matrixGraph, reducer_, blockMatrix_));
    broker_ = brokerptr;
  }

  return(0);
}

//----------------------------------------------------------------------------
int
snl_fei::Factory::createBroker_FEData(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
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


#endif

