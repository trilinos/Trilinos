/*
 * MueLu_AmalgamationInfo_decl.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DECL_HPP_
#define MUELU_AMALGAMATIONINFO_DECL_HPP_

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

namespace MueLu {

/*!
  @class AmalgamationInfo
  @brief minimal container class for storing amalgamation information

  stores maps
   globalamalblockid2myrowid_
   globalamalblockid2globalrowid_

  that are used for unamalgamation
*/

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AmalgamationInfo
    : public BaseClass {
#undef MUELU_AMALGAMATIONINFO_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    AmalgamationInfo() {
      //setObjectLabel(objectLabel);
      globalamalblockid2myrowid_ = Teuchos::null;
      globalamalblockid2globalrowid_ = Teuchos::null;
    }

    virtual ~AmalgamationInfo() {}

    /// Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    void SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid,RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid) const;

    RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > GetMyAmalgamationParams() const;
    RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > GetGlobalAmalgamationParams() const;

  private:

    //! @name amalgamation information variables
    //@{

    /// map: global block id of amalagamated matrix -> vector of local row ids of unamalgamated matrix (only for global block ids of current proc)
    mutable RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid_;   //< used by aggregation factory
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid_; //< used for building overlapping ImportDofMap

    //@}

  };

} // namespace MueLu

#define MUELU_AMALGAMATIONINFO_SHORT
#endif /* MUELU_AMALGAMATIONINFO_DECL_HPP_ */
