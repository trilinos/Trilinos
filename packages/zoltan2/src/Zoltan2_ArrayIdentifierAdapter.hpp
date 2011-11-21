// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_ArrayIdentifierAdapter.hpp

    \brief An input adapter for a simple array of identifiers and weights.
*/


#ifndef _ZOLTAN2_ARRAYIDENTIFIERADAPTER_HPP_
#define _ZOLTAN2_ARRAYIDENTIFIERADAPTER_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierTraits.hpp>

namespace Zoltan2 {

/*! \brief ArrayIdentifierAdapter represents an array of Ids.
 */

template <typename User>
class ArrayIdentifierAdapter : IdentifierAdapter<User> {
private:
public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  ArrayIdentiferAdapter(lno_t numIds, int numWeights, 
    gid_t *ids, scalar_t **weights)
  {
    numIds_ = numIds;
    numWeights_ = numWeights;
    idList_ = ids;
    idWeights_ = weights;
  }

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  std::string inputAdapterName() const {return std::string("ArrayIdentifier");}

  bool haveLocalIds() const {return true;}

  bool haveConsecutiveLocalIds(size_t &base) const 
  {
    base = 0;
    return true;
  }

  ////////////////////////////////////////////////////////////////
  // The IdentifierAdapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIds() const { return numIds_;}

  size_t getIdList(const gid_t *&gids, const lid_t *&lids) const
  {
    gids = idList_;
    lids = NULL;   // it's implied to be 0 to n-1
  }
   
  int getNumWeights() const { return numWeights_; }

  void getWeights(int n, scalar_t *&weights) const
  {
    if ((n>=0) && (n < numWeights))
      weights = idWeights_[n];
    else
      throw std::logic_error("invalid number in getWeights");
  }
};
  
  
}  //namespace Zoltan2
  
#endif
