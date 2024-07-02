// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
