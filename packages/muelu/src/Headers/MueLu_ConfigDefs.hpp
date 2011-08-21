#ifndef MUELU_CONFIGDEFS_HPP
#define MUELU_CONFIGDEFS_HPP

#include "MueLu_config.hpp"

// mem management
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>

//! Namespace for MueLu classes and methods
namespace MueLu {
  // import Teuchos memory management classes into MueLu
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arcp;
  using Teuchos::arcpFromArrayView;
  using Teuchos::rcpFromRef;
  using Teuchos::null;
  using Teuchos::arcp_reinterpret_cast;
  using Teuchos::Array;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcpFromRef;
}

// //! Namespace for MueLu example classes
// namespace MueLuExamples {
// }

#endif /* MUELU_CONFIGDEFS_H */
