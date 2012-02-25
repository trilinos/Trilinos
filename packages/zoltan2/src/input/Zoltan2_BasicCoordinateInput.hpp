// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_BasicCoordinateInput.hpp
    An input adapter for a geometric coordinates (and optional weights)
      that are supplied by the user as pointers to strided arrays.
*/


#ifndef _ZOLTAN2_BASICCOORDINATEINPUT_HPP_
#define _ZOLTAN2_BASICCOORDINATEINPUT_HPP_

#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_StridedInput.hpp>

namespace Zoltan2 {

  /*!  \brief BasicCoordinateInput represents geometric coordinates that are
                supplied by the user as pointers to strided arrays.

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights and coordinates
    \li \c lno_t is the integral data type used by Zoltan2 for local indices and local counts.
    \li \c gno_t is the integral data type used by Zoltan2 to represent global indices and global counts.
    \li \c gid_t is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c gid_t and \c gno_t are the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel, which is used to optimize performance on many-core and multi-core architectures.  If you don't use Kokkos, you can ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent coordinates, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.

  TODO: Global identifiers should be optional.  If the user gives us
    gids in the input adapter, we will include them in the solution.
*/

template <typename User>
  class BasicCoordinateInput : public CoordinateInput<User> {

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef CoordinateInput<User>   base_adapter_t;
  typedef User user_t;

  /*! \brief Constructor
   *
   *  \param dim is the dimension, typically one, two or three, of the
   *                   geometric coordinates.
   *  \param numIds   the local number of coordinates.
   *  \param ids     is a pointer to the coordinate global Ids.  TODO: make
   *                    this optional - if they provide it here, we give it
   *                    back in the solution
   *  \param values a list of \c dim pointers to the coordinate values
   *          corresponding to the \c numIds ids
   *  \param valueStrides a list of \c dim strides for the \c values 
   *        arrays. The k'th coordinate of dimension n can be found at
   *                values[n][k*valueStrides[n]].  If valueStrides
   *              is NULL, it is assumed all strides are one.
   *  \param numWeights the number of weights per coordinate , which may be zero
   *                or greater
   *  \param weights  \c numWeights pointers to arrays of weights.  \c weights
             may be NULL if there are no arrays of weights.
   *  \param weightStrides  a list of \c numWeights strides for the \c weights
   *        arrays. The n'th weight for coordinate k is to be found
   *               at weights[n][k*weightStrides[n]].  If weightStrides
   *              is NULL, it is assumed all strides are one.
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  BasicCoordinateInput( int dim, lno_t numIds, const gid_t *ids, 
    const scalar_t * const *values, int *valueStrides,
    int numWeights, const scalar_t * const * weights, const int *weightStrides):
      env_(rcp(new Environment)), 
      numIds_(numIds), globalNumIds_(), idList_(ids), 
      dimension_(dim), coords_(dim), 
      numWeights_(numWeights), weights_(numWeights)
  {
    typedef StridedInput<lno_t,scalar_t> input_t;

    gno_t tmp = numIds;
    try{
      reduceAll<int, gno_t>(*(env_->comm_), Teuchos::REDUCE_SUM, 1, 
         &tmp, &globalNumIds_);
    }
    Z2_THROW_OUTSIDE_ERROR(*env_, e);

    if (numIds){
      int stride = 1;
      for (int x=0; x < dim; x++){
        if (valueStrides)
          stride = valueStrides[x];
        coords_[x] = rcp<input_t>(new input_t(env_,
          ArrayView<const scalar_t>(values[x], stride*numIds), stride));
      }

      if (numWeights){
        stride = 1;
        for (int w=0; w < numWeights; w++){
          if (weightStrides)
            stride = weightStrides[w];
          weights_[w] = rcp<input_t>(new input_t(env_,
            ArrayView<const scalar_t>(weights[w], stride*numIds), stride));
        }
      }
    }
  }

  /*! Destructor
   */
  ~BasicCoordinateInput() {};

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  std::string inputAdapterName() const {return std::string("BasicCoordinate");}

  size_t getLocalNumberOfObjects() const { return numIds_;}

  int getNumberOfWeightsPerObject() const { return numWeights_;}

  ////////////////////////////////////////////////////
  // The CoordinateInput interface.
  ////////////////////////////////////////////////////

  int getCoordinateDimension() const { return dimension_;}

  int getNumberOfWeights() const { return numWeights_;}  

  size_t getLocalNumberOfCoordinates() const { return numIds_; }

  size_t getGlobalNumberOfCoordinates() const { return globalNumIds_;}

  // TODO make global IDs optional

  size_t getCoordinates(int dim, const gid_t *&gids, const scalar_t *&coords, 
    int &stride) const
  {
    Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid dimension",
      dim >= 0 && dim < dimension_, BASIC_ASSERTION);

    gids = idList_;
    
    size_t length;

    coords_[dim]->getStridedList(length, coords, stride);

    return length;
  }

  size_t getCoordinateWeights(int dim, const scalar_t *&weights, 
    int &stride) const
  {
    Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid dimension",
      dim >= 0 && dim < numWeights_, BASIC_ASSERTION);
    
    size_t length;

    weights_[dim]->getStridedList(length, weights, stride);

    return length;
  }

private:

  // A default environment.  An Environment is an internal Zoltan2
  // class, so input adapters don't usually have one.  But we create
  // one here so we can use it for error handling.

  RCP<const Environment> env_;

  lno_t numIds_;
  gno_t globalNumIds_;

  const gid_t *idList_;

  int dimension_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > coords_;

  int numWeights_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > weights_;
};
  
  
}  //namespace Zoltan2
  
#endif
