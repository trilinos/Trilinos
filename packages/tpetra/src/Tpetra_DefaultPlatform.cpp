#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DefaultPlatform.hpp"

Teuchos::RCP<Tpetra::DefaultPlatform::DefaultPlatformType> Tpetra::DefaultPlatform::platform_ = Teuchos::null;

namespace Tpetra {

  DefaultPlatform::DefaultPlatformType & DefaultPlatform::getDefaultPlatform() {
    if (platform_ == null) {
#ifdef HAVE_TPETRA_MPI
      platform_ = rcp(new MpiPlatform<Kokkos::DefaultNode::DefaultNodeType>());
#else
      platform_ = rcp(new SerialPlatform<Kokkos::DefaultNode::DefaultNodeType>());
#endif
    }
    return *platform_;
  }

}
