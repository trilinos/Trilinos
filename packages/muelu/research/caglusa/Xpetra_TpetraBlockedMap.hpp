#ifndef XPETRA_TPETRABLOCKEDMAP_HPP
#define XPETRA_TPETRABLOCKEDMAP_HPP

#include <Teuchos_RCP.hpp>
#include <Tpetra_BlockedMap_decl.hpp>
#include <Tpetra_BlockedMap_def.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_TpetraVector.hpp>

namespace Xpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class TpetraBlockedMap {
 public:
  using map_type                = Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using lo_vec_type             = Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
  using tpetra_blocked_map_type = Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;

  TpetraBlockedMap(const Teuchos::RCP<const map_type>& pointMap,
                   const Teuchos::RCP<lo_vec_type>& blockSizes) {
    using TpLOVec = TpetraVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
    tpBlockedMap_ = Teuchos::rcp(new tpetra_blocked_map_type(toTpetra(pointMap),
                                                             Teuchos::rcp_dynamic_cast<TpLOVec>(blockSizes)->getTpetra_Vector()));
  }

  RCP<tpetra_blocked_map_type> getTpetra_BlockedMap() const {
    return tpBlockedMap_;
  }

 private:
  RCP<tpetra_blocked_map_type> tpBlockedMap_;
};

}  // namespace Xpetra

#endif  // XPETRA_TPETRABLOCKEDMAP_HPP
