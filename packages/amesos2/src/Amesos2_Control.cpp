/**
  \file   Amesos2_Control.cpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Sat Jan 16 09:23:38 2010
  
  \brief  Implementation for Amesos2::Control
*/

#ifndef AMESOS2_CONTROL_CPP
#define AMESOS2_CONTROL_CPP

#include <sstream>
#include "Amesos2_Control.hpp"

namespace Amesos {


void Control::setControlParameters(
  const Teuchos::RCP<Teuchos::ParameterList> & parameterList )
{
  // Add this value to all diagonal elements which are structurally
  // non-zero. No change is made to non-zero structure of the matrix.
  if( parameterList->isParameter("AddToDiag") ){
    addToDiag_ = parameterList->get<double>("AddToDiag");
  }
  
  // Add zero to diagonal if diagonal element is not present.
  // - not supported for matrices which are missing elements from the diagonal.  See bug #1928 for discussion
  if( parameterList->isParameter("AddZeroToDiag") ){
    addZeroToDiag_ = parameterList->get<bool>("AddZeroToDiag");
  }

  // Matrix property, defined internally in Amesos2::Mumps as an integer,
  // whose value can be:
  // - 0 : general unsymmetric matrix;
  // - 1 : SPD;
  // - 2 : general symmetric matrix.
  if( parameterList->isParameter("MatrixProperty") ) {
    std::string MatrixProperty;
    MatrixProperty = parameterList->get<std::string>("MatrixProperty");
    if( MatrixProperty == "general" )
      matrixProperty_ = 0;
    else if( MatrixProperty == "SPD" )   
      matrixProperty_ = 1;
    else if( MatrixProperty == "symmetric" ) 
      matrixProperty_ = 2;
    else {
      std::ostringstream oss;
      oss << "Amesos2 : ERROR" << std::endl 
          << "Amesos2 : matrixProperty value not recognized ("
          << MatrixProperty << ")" << std::endl;
      TEST_FOR_EXCEPTION(true,
        std::runtime_error,
        oss.str());
    }

    // Threshold for determining if refactorize worked OK
    // UNUSED at present - KSS June 2004
    if( parameterList->isParameter("RcondThreshold") ){
      rcond_threshold_ = parameterList->get<double>("RcondThreshold");
    }

    // Determine whether to Refactorize.
    if( parameterList->isParameter("Refactorize") ){
      refactorize_ = parameterList->get<bool>("Refactorize");
    }

    // Define how many processes to use in the ScaLAPACK factor and
    // solve.  If (-1), a heuristic is used to determine the number of
    // processes to use
    if( parameterList->isParameter("MaxProcs") ){
      maxProcesses_ = parameterList->get<int>("MaxProcs");
    }
    
    // Scaling method: 0: none, 1: use method's default, 2: use
    // the method's 1st alternative, 3: etc.
    if( parameterList->isParameter("ScaleMethod") ){
      scaleMethod_ = parameterList->get<int>("ScaleMethod");
    }

    // If true, the Amesos2 class should reindex the matrix to
    // standard indexing (i.e. 0-(n-1)).  At present, only
    // Amesos2::Klu supports this option.
    if( parameterList->isParameter("Reindex") ){
      reindex_ = parameterList->get<bool>("Reindex");
    }
  }
}


} // end namespace Amesos

#endif	// AMESOS2_CONTROL_CPP
