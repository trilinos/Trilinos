#include "Cthulhu_ConfigDefs.hpp"

#ifdef HAVE_CTHULHU_EPETRA

#include "Cthulhu_EpetraMap.hpp"

namespace Cthulhu {

  const Epetra_Map & toEpetra(const Cthulhu::Map<int,int> &map) {
    const Cthulhu::EpetraMap & epetraMap = dynamic_cast<const Cthulhu::EpetraMap &>(map);
    return epetraMap.getEpetra_Map();
  }

}

// /** \brief  Returns true if \c map is identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
//     \relates Cthulhu::EpetraMap */
// bool operator== (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) {
//  return map1.isSameAs(map2);
// }

// /** \brief Returns true if \c map is not identical to this map. Implemented in Cthulhu::EpetraMap::isSameAs().
//     \relates Cthulhu::EpetraMap */
// bool operator!= (const Cthulhu::EpetraMap &map1, const Cthulhu::EpetraMap &map2) {
//  return !map1.isSameAs(map2);
// }

#endif
