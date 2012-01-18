// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_CoordinateInput.hpp

    \brief The abstract interface for an input adapter representing geometric
               coordinates with optional weights.
*/

#ifndef _ZOLTAN2_COORDINATEINPUT_HPP_
#define _ZOLTAN2_COORDINATEINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*!  \brief CoordinateInput defines the interface for input of geometric 
                coordinates with optional weights.

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights and coordinates.
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
*/

template <typename User>
  class CoordinateInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  /*! \brief Pure virtual destructor
   */
  virtual ~CoordinateInput() {};

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  enum InputAdapterType inputAdapterType() const {return CoordinateAdapterType;}

  ////////////////////////////////////////////////////
  // My interface.
  ////////////////////////////////////////////////////

  /*! \brief Return dimension of the coordinates.
   *   \return the number of coordinates (typically one, two or three).
   */
  virtual int getCoordinateDimension() const = 0;


  /*! \brief Return the number of weights per coordinate.
   *   \return the count of weights, zero or more per coordinate.
   */
  virtual int getNumberOfWeights() const = 0;

  /*! Return the number of coordinates on this process.
   *   \return  the count of coordinates on the local process.
   */
  virtual size_t getLocalNumberOfCoordinates() const = 0;

  /*! Return the number of coordinates in the entire problem.
   *   \return  the global count of coordinates.
   */
  virtual size_t getGlobalNumberOfCoordinates() const = 0;

  /*! Provide a pointer to one dimension of this process' coordinates.
      \param dim  is a value from 0 to one less than 
         getLocalNumberOfCoordinates() specifying which dimension is
         being provided in the coords list.
      \param coords  points to a list of coordinate values for the dimension.
      \param stride  describes the layout of the coordinate values in
              the coords list.  If stride is one, then the ith coordinate
              value is coords[i], but if stride is two, then the
              ith coordinate value is coords[2*i].

       \return The length of the \c coords list.  This may be more than
              getLocalNumberOfCoordinates() because the \c stride
              may be more than one.

      TODO make global IDs optional - we'll return them
        in the solution if they include then in this call.

      Zoltan2 does not copy your data.  The data pointed to coords
      must remain valid for the lifetime of this InputAdapter.
   */

  virtual size_t getCoordinates(int dim, const gid_t *&gids, 
    const scalar_t *&coords, int &stride) const = 0;

  /*! \brief  Provide a pointer to the weights, if any, corresponding 
       to the coordinates returned in getCoordinates(). 

      \param dimension ranges from zero to one less than getNumberOfWeights()
      \param weights is the list of weights of the given dimension for
           the coordinates listed in getCoordinates().
       \param stride The k'th weight is located at weights[stride*k]
       \return The number of weights listed, which should be the same
                  as the number of elements listed in getCoordinates().
   */

  virtual size_t getCoordinateWeights(int dimension,
     const scalar_t *&weights, int &stride) const = 0;

  /*! \brief Apply a PartitioningSolution to an input.
   *
   *  This is not a required part of the CoordinateInput interface. However
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
