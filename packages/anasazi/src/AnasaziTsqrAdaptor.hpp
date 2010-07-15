#ifndef __Anasazi_TsqrAdaptor_hpp
#define __Anasazi_TsqrAdaptor_hpp

#include "TsqrAdaptor.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Anasazi {
  template< class ScalarType, class MultiVectorType >
  class TsqrAdaptor< class ScalarType, class MultiVectorType > 
  {
  public:
    typedef TsqrAdaptor< typename MultiVectorType::scalar_type,
			 typename MultiVectorType::local_ordinal_type,
			 typename MultiVectorType::global_ordinal_type,
			 typename MultiVectorType::node_type > adaptor_type;
  };

  // FIXME mfh 14 Jul 2010: this belongs in ../tpetra/src
  template< class S, class LO, class GO, class Node >
  class TsqrAdaptor< ScalarType, Tpetra::MultiVector< S, LO, GO, Node > >
  {
  public:
    typedef TsqrTpetraAdaptor< S, LO, GO, Node > adaptor_type;
  };

  // FIXME mfh 14 Jul 2010: this belongs in ../epetra/src
  template <>
  class TsqrAdaptor< double, Epetra_MultiVector >
  {
  public:
    typedef TsqrEpetraAdaptor adaptor_type;
  };

} // namespace Anasazi

#endif // __Anasazi_TsqrAdaptor_hpp
