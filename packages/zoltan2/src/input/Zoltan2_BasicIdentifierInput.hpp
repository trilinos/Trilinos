// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_BasicIdentifierInput.hpp
    \brief Defines the BasicIdentifierInput class.
*/

#ifndef _ZOLTAN2_BASICIDENTIFIERINPUT_HPP_
#define _ZOLTAN2_BASICIDENTIFIERINPUT_HPP_

#include <Zoltan2_IdentifierInput.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*! \brief This class represents a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  The user supplies the identifiers and weights by way of pointers
 *    to arrays.  
 */

template <typename User>
class BasicIdentifierInput: public IdentifierInput<User> {

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef IdentifierInput<User>       base_adapter_t;
  typedef User user_t;

  /*! \brief Constructor
   *  \param numIds is the number of identifiers in the list
   *  \param ids should point to a list of numIds identifiers.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per identifier is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  BasicIdentifierInput( lno_t numIds, const gid_t *idPtr, 
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  string inputAdapterName() const {return string("BasicIdentifier");}

  size_t getLocalNumberOfObjects() const { return numIds_;}

  int getNumberOfWeightsPerObject() const { return weights_.size();}

  ////////////////////////////////////////////////////////////////
  // The IdentifierInput interface.
  // This is the interface that would be called by a model or a problem .
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumberOfIdentifiers() const { return numIds_;}
   
  int getNumberOfWeights() const { return weights_.size(); }

  size_t getIdentifierList(const gid_t *&Ids) const
  {
    Ids = idList_;
    return numIds_;
  }

  size_t getIdentifierWeights(int dimension,
     const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid weight dimension",
      dimension >= 0 && dimension < weights_.size(), BASIC_ASSERTION);

    size_t length;
    weights_[dimension]->getStridedList(length, weights, stride);
    return length;
  }

private:

  RCP<const Environment> env_;
  lno_t numIds_;
  const gid_t *idList_;
  Array<RCP<StridedData<lno_t, scalar_t> > > weights_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
  BasicIdentifierInput<User>::BasicIdentifierInput(
    lno_t numIds, const gid_t *idPtr,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      numIds_(numIds), idList_(idPtr), weights_(weights.size())
{
  env_ = rcp(new Environment);    // for error messages
  size_t numWeights = weights.size();

  if (numWeights){
    typedef StridedData<lno_t,scalar_t> input_t;
    if (weightStrides.size())
      for (int i=0; i < numWeights; i++)
        weights_[i] = rcp<input_t>(new input_t(
          ArrayView<const scalar_t>(weights[i], weightStrides[i]*numIds),
          weightStrides[i]));
    else
      for (int i=0; i < numWeights; i++)
        weights_[i] = rcp<input_t>(new input_t(
          ArrayView<const scalar_t>(weights[i], numIds), 1));
  }
}

  
  
}  //namespace Zoltan2
  
#endif
