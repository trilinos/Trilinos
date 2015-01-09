// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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

#ifndef __Anasazi_TsqrAdaptor_hpp
#define __Anasazi_TsqrAdaptor_hpp

#include "TsqrAdaptor.hpp"
#include "TsqrRandomizer.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Anasazi {
  /// \class TsqrAdaptor
  /// \brief Map from multivector class to TSQR adaptor class
  template< class ScalarType, class MultiVectorType >
  class TsqrAdaptor
  {
  public:
    typedef TSQR::Trilinos::TsqrAdaptor< typename MultiVectorType::scalar_type,
					 typename MultiVectorType::local_ordinal_type,
					 typename MultiVectorType::global_ordinal_type,
					 typename MultiVectorType::node_type > adaptor_type;
  };

  // FIXME mfh 14 Jul 2010: this belongs in ../epetra/src
  // template <>
  // class TsqrAdaptor< double, Epetra_MultiVector >
  // {
  // public:
  //   typedef TsqrEpetraAdaptor adaptor_type;
  // };

} // namespace Anasazi

#endif // __Anasazi_TsqrAdaptor_hpp
