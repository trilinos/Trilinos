// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_VectorInput.hpp

    \brief The interface for an input adapter for a 
    distributed vector with optional weights. 
*/


#ifndef _ZOLTAN2_VECTORINPUT_HPP_
#define _ZOLTAN2_VECTORINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2 {

  /*!  \brief VectorInput defines the interface for vector input adapters.

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights and vector element values.
    \li \c lno_t is the integral data type used by Zoltan2 for local indices and local counts.
    \li \c gno_t is the integral data type used by Zoltan2 to represent global indices and global counts.
    \li \c gid_t is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c gid_t and \c gno_t are the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel, which is used to optimize performance on many-core and multi-core architectures.  If you don't use Kokkos, you can ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent a vector, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.

    VectorInput may be a single vector or a set of corresponding vectors
    which have with the
    same global identifiers and the same distribution across processes.
    (For example, there is a Trilinos Xpetra::Multivector input adapter
    which is a sub class of VectorInput.)
 
*/

template <typename User>
  class VectorInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  /*! Pure virtual destructor
   */
  virtual ~VectorInput() {};

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  enum InputAdapterType inputAdapterType() const {return VectorAdapterType;}

  ////////////////////////////////////////////////////
  // My interface.
  ////////////////////////////////////////////////////

  /*! \brief Return the number of vectors (typically one).
   */
  virtual int getNumberOfVectors() const = 0;

  /*! \brief Return the number of weights per vector element.
   */
  virtual int getNumberOfWeights() const = 0;

  /*! \brief Return the length of the portion of the vector on this process.
   */
  virtual size_t getLocalLength() const = 0;

  /*! \brief Return the global length of the vector on this process.
   */
  virtual size_t getGlobalLength() const = 0;

  /*! \brief Provide a pointer to the vertex elements.  If the VectorInput
       represents more than one vector, vector zero is implied.

      \param ids will on return point to the list of global Ids for 
        each element on this process.  TODO: make Ids optional.
      \param elements will on return point to the vector values
        corresponding to the global Ids.
      \param stride the k'th element is located at elements[stride*k]
      \return The number of ids in the Ids list.
   */

  virtual size_t getVector(const gid_t *&ids, 
     const scalar_t *&elements, int &stride) const = 0;

  /*! \brief Provide a pointer to the elements of the specified vector.

      \param i ranges from zero to one less than getNumberOfVector(), and
         represents the vector for which data is being requested.
      \param ids will on return point to the list of global Ids for 
        each element on this process.  TODO: make Ids optional.
      \param elements will on return point to the vector values
        corresponding to the global Ids.
      \param stride the k'th element is located at elements[stride*k]
      \return The number of ids in the Ids list.
   */

  virtual size_t getVector(int i, const gid_t *&ids, 
     const scalar_t *&elements, int &stride) const = 0;

  /*! \brief  Provide a pointer to the weights corresponding to the elements
       returned in getVector().
        
      \param dimension ranges from zero to one less than getNumberOfWeights()
      \param weights is the list of weights of the given dimension for
           the elements listed in getVector.
       \param stride The k'th weight is located at weights[stride*k]
       \return The number of weights listed, which should be the same
                  as the number of elements listed in getVector().
   */

  virtual size_t getVectorWeights(int dimension,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the VectorInput interface. However
   *  if the Caller calls a Problem method to redistribute data, it needs
   *  this method to perform the redistribution.
   *
   *  \param in  An input object with a structure and assignment of
   *           of global Ids to processes that matches that of the input
   *           data that instantiated this InputAdapter.
   *  \param out On return this should point to a newly created object
   *            with the specified partitioning.
   *  \param solution  The Solution object created by a Problem should
   *      be supplied as the third argument.  It must have been templated
   *      on user data that has the same global ID distribution as this
   *      user data.
   *  \return   Returns the number of local Ids in the new partitioning.
   */

  template <typename User2>
    size_t applyPartitioningSolution(User &in, User *&out,
         const PartitioningSolution<User2> &solution)
  {
    return 0;
  } 
};
  
  
}  //namespace Zoltan2
  
#endif
