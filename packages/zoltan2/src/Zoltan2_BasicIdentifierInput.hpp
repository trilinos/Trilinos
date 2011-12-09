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
#include <Zoltan2_StridedInput.hpp>

namespace Zoltan2 {

/*! \brief This class represents a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  A pointer to the global identifiers is supplied in the constructor.
 */

template <typename User>
class BasicIdentifierInput: IdentifierInput<User> {

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef IdentifierInput<User>       base_adapter_t;

  /*! Constructor
      \param numIds is the number of identifiers in the list
      \param numWeights is the number of weights provided for each
                        identifier.  Weights are optional.
      \param ids should point to a list of numIds identifiers.
      \param weights should point to a list of numWeights pointers.  Each
          pointer should point to a list of weights for the identifiers.
      \param strides should point to a list of numWeights integers. 
          strides[i] is the stride for list weights[i].  If strides is
          NULL, it will be assumed that each stride is one.
   */

  BasicIdentifierInput(lno_t numIds, int numWeights, const gid_t *ids, 
    const scalar_t **weights, const int *strides): 
      numIds_(numIds), idList_(ids), weights_(numWeights)
  {
    if (numWeights){
      RCP<const Environment> env = rcp(new Environment); 
      if (strides)
        for (int i=0; i < numWeights; i++)
          weights_[0] = StridedInput(env, weights[i], strides[i]);
      else
        for (int i=0; i < numWeights; i++)
          weights_[0] = StridedInput(env, weights[i], 1);
    }
  }

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  std::string inputAdapterName() const {return std::string("BasicIdentifier");}

  ////////////////////////////////////////////////////////////////
  // The IdentifierInput interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIds() const { return numIds_;}
   
  int getNumWeights() const { return weights_.size(); }

  size_t getIdList(const gid_t *&Ids, const scalar_t **&weights, 
    const int *&strides) const
  {
    gids = idList_;

    int len, stride;
    const scalar_t *vec;

    for (int i=0; i < numWeights; i++){
      weights_[i].getStridedList(len, vec, stride);
      weights[i] = vec
      strides[i] = stride;
    }

    return numIds_;
  }



private:

  lno_t numIds_;
  const gid_t *idList_;
  std::vector<StridedInput<scalar_t> > weights_;
};
  
  
}  //namespace Zoltan2
  
#endif
