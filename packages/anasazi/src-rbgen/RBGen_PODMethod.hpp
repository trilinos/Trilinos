// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_POD_METHOD_HPP
#define RBGEN_POD_METHOD_HPP

#include "RBGen_ConfigDefs.h" 

namespace RBGen {

  //! Abstract base class for reduced basis POD methods
  template<class ScalarType>  
  class PODMethod {
    
  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    PODMethod() {};

    //! Destructor.
    virtual ~PODMethod() {};
    //@}

    //! @name Get methods
    //@{
    
    //! Returns the singular values computed corresponding to the reduced basis.
    virtual std::vector<ScalarType> getSingularValues() const = 0;

    //@}
    
  };
  
} // end of RBGen namespace

#endif // RBGEN_POD_METHOD_HPP
