// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierInput.hpp

    \brief The abstract interface for an input adapter that simply
             represents a list of identifiers with optional weights.
*/


#ifndef _ZOLTAN2_IDENTIFIERINPUT_HPP_
#define _ZOLTAN2_IDENTIFIERINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*!  \brief IdentifierInput defines the interface for input adapters 
           that represent a list of identifiers and weights.

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights .
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
*/

template <typename User>
  class IdentifierInput : public InputAdapter {

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  /*! Pure virtual destructor
   */
  virtual ~IdentifierInput() {};

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  enum InputAdapterType inputAdapterType() const {return IdentifierAdapterType;}

  ////////////////////////////////////////////////////
  // My interface.
  ////////////////////////////////////////////////////

  /*! \brief Return the number of identifiers on this process.
   */
  virtual size_t getLocalNumberOfIdentifiers() const = 0;

  /*! \brief Return the number of weights associated with each identifier.
   */
  virtual int getNumberOfWeights() const = 0;

  /*! \brief Provide pointers to this process' identifiers and optional weights.

      \param Ids will on return point to the list of the global Ids for 
        this process.
      \param weights on return will contain numWeights pointers, each
         pointing to one of the lists of weights.  (The caller must
         allocate the array of numWeighst pointers.)
      \param strides on return will contain numWeights numbers, one for each
         weight list, indicating the stride of that list.  
         (The caller must allocate the array of numWeights integers.)

       \return The number of ids in the Ids list.
   */

  virtual size_t getIdentifierList(gid_t const **Ids, scalar_t const **weights,
    int *strides) const = 0;


 /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the IdentifierInput interface. However
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
