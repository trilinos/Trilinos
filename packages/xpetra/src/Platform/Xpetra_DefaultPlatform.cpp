#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"

Teuchos::RCP<Xpetra::DefaultPlatform::DefaultPlatformType> Xpetra::DefaultPlatform::platform_ = Teuchos::null;

namespace Xpetra {

  DefaultPlatform::DefaultPlatformType &DefaultPlatform::getDefaultPlatform() {
    if (!platform_.get()) {
#ifdef HAVE_MPI
      platform_ = Teuchos::rcp(new MpiPlatform<Kokkos::DefaultNode::DefaultNodeType>(Kokkos::DefaultNode::getDefaultNode()));
#else
      platform_ = Teuchos::rcp(new SerialPlatform<Kokkos::DefaultNode::DefaultNodeType>(Kokkos::DefaultNode::getDefaultNode()));
#endif
    }
    return *platform_;
  }

}
