#ifndef AMESOS2_MATRIXTRAITS_HPP
#define AMESOS2_MATRIXTRAITS_HPP

#include "Amesos2_config.h"

#include <Tpetra_CrsMatrix.hpp>

#ifdef HAVE_AMESOS2_EPETRA
#  include <Tpetra_DefaultPlatform.hpp>
#  include <Epetra_RowMatrix.h>
#  include <Epetra_CrsMatrix.h>
// #  include <Epetra_MsrMatrix.h>
#  include <Epetra_VbrMatrix.h>
// and perhaps some others later...
#endif

#include "Amesos2_Util.hpp"

namespace Amesos {

  // The declaration
  template <class Matrix>
  struct MatrixTraits {};

  /*******************
   * Specializations *
   *******************/

  template < typename Scalar,
	     typename LocalOrdinal,
	     typename GlobalOrdinal,
	     typename Node >
  struct MatrixTraits<
    Tpetra::RowMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef GlobalOrdinal global_ordinal_t;
    typedef Node node_t;

    typedef Util::row_access major_access;
  };

  template < typename Scalar,
	     typename LocalOrdinal,
	     typename GlobalOrdinal,
	     typename Node,
	     typename LocalMatOps >
  struct MatrixTraits<
    Tpetra::CrsMatrix<Scalar,
		      LocalOrdinal,
		      GlobalOrdinal,
		      Node,
		      LocalMatOps> > {
    typedef Scalar scalar_t;
    typedef LocalOrdinal local_ordinal_t;
    typedef GlobalOrdinal global_ordinal_t;
    typedef Node node_t;

    typedef Util::row_access major_access;
  };

#ifdef HAVE_AMESOS2_EPETRA

  template <>
  struct MatrixTraits<Epetra_RowMatrix> {
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef int global_ordinal_t;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_t;
    
    typedef Util::row_access major_access;
  };

  template <>
  struct MatrixTraits<Epetra_CrsMatrix> {
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef int global_ordinal_t;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_t;
    
    typedef Util::row_access major_access;
  };

  // template <>
  // struct MatrixTraits<Epetra_MsrMatrix> {
  //   typedef double scalar_t;
  //   typedef int local_ordinal_t;
  //   typedef int global_ordinal_t;
  //   typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_t;
    
  //   typedef Util::row_access major_access;
  // };

  template <>
  struct MatrixTraits<Epetra_VbrMatrix> {
    typedef double scalar_t;
    typedef int local_ordinal_t;
    typedef int global_ordinal_t;
    typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_t;
    
    typedef Util::row_access major_access;
  };

#endif

};

#endif	// AMESOS2_MATRIXTRAITS_HPP
