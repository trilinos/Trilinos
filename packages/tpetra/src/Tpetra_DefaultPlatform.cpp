#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"

Teuchos::RCP<Tpetra::DefaultPlatform::DefaultPlatformType> Tpetra::DefaultPlatform::platform_ = Teuchos::null;

namespace Tpetra {

  DefaultPlatform::DefaultPlatformType &DefaultPlatform::getDefaultPlatform() {
    if (!platform_.get()) {
#ifdef HAVE_TPETRA_MPI
      platform_ = Teuchos::rcp(new MpiPlatform<Kokkos::DefaultNode::DefaultNodeType>(Kokkos::DefaultNode::getDefaultNode()));
#else
      platform_ = Teuchos::rcp(new SerialPlatform<Kokkos::DefaultNode::DefaultNodeType>(Kokkos::DefaultNode::getDefaultNode()));
#endif
    }
    return *platform_;
  }

}
