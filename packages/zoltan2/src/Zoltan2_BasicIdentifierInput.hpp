// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_BasicIdentifierInput.hpp

    \brief An input adapter for a simple array of identifiers and weights.
*/


#ifndef _ZOLTAN2_BASICIDENTIFIERINPUT_HPP_
#define _ZOLTAN2_BASICIDENTIFIERINPUT_HPP_

#include <Zoltan2_IdentifierInput.hpp>

namespace Zoltan2 {

/*! \brief This class represents a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  A pointer to the global identifiers is supplied in the constructor.
 *  Local identifiers are implied to be consecutive beginning at zero.
 */

template <typename User>
class BasicIdentifierInput: IdentifierInput<User> {

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef IdentifierInput<User>       base_adapter_t;

  BasicIdentifierInput(lno_t numIds, int numWeights, const gid_t *ids, 
    const scalar_t *weights): numIds_(numIds), numWeights_(numWeights), 
      idList_(ids), idWeights_(weights) { }

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  std::string inputAdapterName() const {return std::string("BasicIdentifier");}

  bool haveLocalIds() const {return true;}

  bool haveConsecutiveLocalIds(size_t &base) const 
  {
    base = 0;
    return true;
  }

  ////////////////////////////////////////////////////////////////
  // The IdentifierInput interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIds() const { return numIds_;}

  size_t getIdList(const gid_t *&gids, const lid_t *&lids,
    const scalar_t *&weights) const
  {
    gids = idList_;
    lids = NULL;   // it's implied to be 0 to n-1
    weights = idWeights_;
    return numIds_;
  }
   
  int getNumWeights() const { return numWeights_; }

private:

  lno_t numIds_;
  int numWeights_;
  const gid_t *idList_;
  const scalar_t *idWeights_;
};
  
  
}  //namespace Zoltan2
  
#endif
