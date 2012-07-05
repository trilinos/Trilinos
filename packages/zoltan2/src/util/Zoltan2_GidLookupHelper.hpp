// @HEADER
// ***********************************************************************
//            copyright
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GidLookupHelper.hpp
   \brief Defines GidLookupHelper class.
*/

#ifndef _ZOLTAN2_GIDLOOKUPHELPER
#define _ZOLTAN2_GIDLOOKUPHELPER

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_IdentifierTraits.hpp>
#include <Teuchos_Hashtable.hpp>
#include <map>

namespace Zoltan2
{
/*! \brief A class to do GID lookups.
 *
 * Processing user global IDs may require GID -> index look ups.
 * Some T support generation of unique double hash keys for any value.
 * Ordinals are an example.  So do not, such as strings.
 *
 * If double hash keys are supported, a hash table will be created.  If
 * not, a std::map will be created.
 */
template <typename T, typename lno_t>
  class GidLookupHelper{

private:

  RCP<const Environment> env_;
  ArrayRCP<const T> gidList_;
  bool useHashTable_;

  map<T, lno_t> indexMap_;

  RCP<Teuchos::Hashtable<double, lno_t> > indexHash_;

public:
  /*! \brief Constructor
   */
  GidLookupHelper(const RCP<const Environment> &env, 
    const ArrayRCP<const T> &gidList);

  /*! \brief Construct an empty lookup helper.
   */
  GidLookupHelper();

  /*! \brief the index of the global id in the \c gidList passed
   *   into the constructor.
   *
   *  If duplicate gids appear in the list, it is the location of
   *  the first of these.
   */
  lno_t lookup(const T gid) const;

  /*! \brief the number of unique gids in the list.
   */
  lno_t size() const {
    if (useHashTable_)
      return indexHash_->size();
    else
      return indexMap_.size();
  }

  /*! \brief  Return a list of the indices in gidList that were
   *              included in the lookup table (all if there were
   *              no duplicates.)
   *  Returns indices from lowest to highest
   */
  void getIndices(ArrayRCP<lno_t> &indices) const;
};

template<typename T, typename lno_t>
  GidLookupHelper<T, lno_t>::GidLookupHelper():
    env_(rcp(new Environment)), gidList_(), useHashTable_(false), 
    indexMap_(), indexHash_()
 {}

template<typename T, typename lno_t>
  GidLookupHelper<T, lno_t>::GidLookupHelper(
    const RCP<const Environment> &env, 
    const ArrayRCP<const T> &gidList):
    env_(env),
    gidList_(gidList), useHashTable_(false), indexMap_(), indexHash_()
{
  lno_t len = gidList_.size();

  if (len < 1)
    return;

  if (IdentifierTraits<T>::hasUniqueKey())
    useHashTable_ = true;

  const T *ids = gidList_.getRawPtr();

  if (!useHashTable_){
    try{
      for (lno_t i=0; i < gidList.size(); i++){
        typename map<T, lno_t>::iterator rec = indexMap_.find(*ids);
        if (rec == indexMap_.end())
          indexMap_[*ids] = i;
        ids++;
      }
    }
    catch (const std::exception &e){
      env_->localMemoryAssertion(__FILE__, __LINE__, len, false);
    }
  }
  else{
    typedef typename Teuchos::Hashtable<double, lno_t> id2index_hash_t;
    id2index_hash_t *p = NULL;
  
    try{
      p = new id2index_hash_t(len);
    }
    catch (const std::exception &e){
      env_->localMemoryAssertion(__FILE__, __LINE__, len, false);
    }
  
    for (lno_t i=0; i < len; i++){
      double key = IdentifierTraits<T>::key(*ids++);
      try{
        if (!p->containsKey(key))
          p->put(key, i);
      }
      Z2_THROW_OUTSIDE_ERROR(*env_);
    }

    indexHash_ = rcp<id2index_hash_t>(p);
  }
}

template<typename T, typename lno_t>
  lno_t GidLookupHelper<T, lno_t>::lookup(const T gid) const
{
  lno_t idx;
  bool badId=false;
  if (useHashTable_){
    try{
      double key = IdentifierTraits<T>::key(gid);
      idx = indexHash_->get(key);
    }
    catch (const std::exception &e) {
      badId = true;
    }
  }
  else{
    typename map<T, lno_t>::const_iterator next = indexMap_.find(gid);
    if (next == indexMap_.end())
      badId = true;
    else
      idx = (*next).second;
  }

  env_->localInputAssertion(__FILE__, __LINE__, "invalid global id",
    badId==false, BASIC_ASSERTION);

  return idx;
}

template<typename T, typename lno_t>
  void GidLookupHelper<T, lno_t>::getIndices(
    ArrayRCP<lno_t> &indices) const
{
  lno_t numIds = size();
  lno_t *idx = new lno_t [numIds];
  ArrayRCP<lno_t> indexList = arcp(idx, 0, numIds, true);

  if (numIds == gidList_.size()){
    for (lno_t i=0; i < gidList_.size(); i++){
       *idx++ = i;
    }
  }
  else{
    const T *ids = gidList_.getRawPtr();
    set<T> allGids;

    for (lno_t i=0; i < gidList_.size(); i++){
      typename set<T>::iterator rec = allGids.find(ids[i]);
      if (rec == allGids.end()){
        allGids.insert(ids[i]);
        *idx++ = i;
      }
    }
  }

  indices = indexList;
}

} // namespace Z2

#endif // _ZOLTAN2_IDENTIFIERTRAITS
