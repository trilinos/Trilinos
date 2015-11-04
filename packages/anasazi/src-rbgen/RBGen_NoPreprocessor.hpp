// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef RBGEN_NO_PREPROCESSOR_HPP
#define RBGEN_NO_PREPROCESSOR_HPP

#include "RBGen_Preprocessor.hpp"

namespace RBGen {

  //! Null preprocessor concrete class.
  template< class DataSetType > 
  class NoPreprocessor : public Preprocessor< DataSetType > {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    NoPreprocessor() {};

    //! Destructor.
    virtual ~NoPreprocessor() {};
    //@}

    //! @name Initialization/Reset Methods
    //@{

    //! Initialize preprocessor
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params, 
                     const Teuchos::RCP< FileIOHandler<DataSetType> >& fileio ) {};

    //! Reset preprocessor
    void Reset() {};
    //@}

    //! @name Preprocess Methods
    //@{

    //! Preprocess the snapshot set passed in
    void Preprocess( Teuchos::RCP<DataSetType>& ss ) {};
    //@}

    //! @name Status Methods
    //@{

    //! Return initialized status of the preprocessor
    bool isInitialized() const { return true; };
    //@}
  };

} // end of RBGen namespace

#endif // RBGEN_NO_PREPROCESSOR_HPP
