#include "MueLu_MutuallyExclusiveTime.hpp"
#include "MueLu_FactoryBase.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  template <class TagName>
  std::stack<MutuallyExclusiveTime<TagName>*> MutuallyExclusiveTime<TagName>::timerStack_;

  //FIXME: move this:
  template class MutuallyExclusiveTime<FactoryBase>;
  template class MutuallyExclusiveTime<Level>;

} // namespace MueLu
