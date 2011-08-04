#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_DefaultPlatform.hpp"

Teuchos::RCP<Cthulhu::DefaultPlatform::DefaultPlatformType> Cthulhu::DefaultPlatform::platform_ = Teuchos::null;

namespace Cthulhu {

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
