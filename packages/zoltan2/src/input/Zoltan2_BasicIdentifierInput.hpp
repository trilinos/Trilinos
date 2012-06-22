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
 *
    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent coordinates, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.

    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.

 */

template <typename User>
  class BasicIdentifierInput: public IdentifierInput<User> {

public:

  typedef typename InputTraits<User>::scalar_t    scalar_t;
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

  size_t getObjectWeights(int dim, const scalar_t *&wgt, int &stride) const
  {
    return getIdentifierWeights(dim, wgt, stride);
  }

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
    weights_[dimension].getStridedList(length, weights, stride);
    return length;
  }

private:

  RCP<const Environment> env_;
  lno_t numIds_;
  const gid_t *idList_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
  BasicIdentifierInput<User>::BasicIdentifierInput(
    lno_t numIds, const gid_t *idPtr,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      numIds_(numIds), idList_(idPtr), weights_()
{
  typedef StridedData<lno_t,scalar_t> input_t;
  env_ = rcp(new Environment);    // for error messages
  size_t numWeights = weights.size();

  if (numWeights > 0){
    weights_ = arcp(new input_t [numWeights], 0, numWeights, true);

    if (numIds > 0){
      for (size_t i=0; i < numWeights; i++){
        int stride = weightStrides.size() ? weightStrides[i] : 1;
        ArrayRCP<const scalar_t> wgtV(weights[i], 0, stride*numIds, false);
        weights_[i] = input_t(wgtV, stride);
      }
    }
  }
}

  
  
}  //namespace Zoltan2
  
#endif
