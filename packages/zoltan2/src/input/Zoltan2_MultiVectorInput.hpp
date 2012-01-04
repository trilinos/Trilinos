// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MultiVectorInput.hpp

    \brief The abstract interface for an input adapter representing a 
    set of distributed vectors with each vector in the set having the
     same mapping from index to process.
*/

#ifndef _ZOLTAN2_MULTIVECTORINPUT_HPP_
#define _ZOLTAN2_MULTIVECTORINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <string>

namespace Zoltan2 {

/*! Zoltan2::MultiVectorInput
    \brief MultiVectorInput defines the interface for input adapters sets of vectors

    GID: the type for the application's global Ids
    LNO: the integral type that Zoltan2 will use for local counters.
    GNO: the integral type that Zoltan2 will use for the global 
      counts and identifiers.  It needs to be large enough for the
      problem's number of objects.

    The template parameter is the the User's vector data structure.
    The user must have previously defined traits for their vector 
    data structure, supplying Zoltan2 mainly with their data types.
 
    There is no reason that we could not use MultiVectorInput
    to represent Vector also.  However this may
    be confusing to people with single vector input who
     would have to derive from MultiVectorInput.
*/

template <typename User>
  class MultiVectorInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  enum InputAdapterType inputAdapterType() const {return MultiVectorAdapterType;}

  /*! Pure virtual destructor
   */
  virtual ~MultiVectorInput() {};

  /*! Return the number of vectors in the multi-vector
   */
  virtual size_t getNumVectors() const = 0;

  /*! Return the length of the portion of the vector on this process.
   */
  virtual size_t getLocalLength() const = 0;

  /*! Return the global length of the vector on this process.
   */
  virtual size_t getGlobalLength() const = 0;

  /** TODO getStride - what exactly does this mean
   */

  /*! Sets pointers to this process' vertex elements.
      \param i  the vector to return, i ranges from 0 to getNumVectors()-1
      \param Ids will on return point to the list of the global Ids for 
        each element on this process.
      \param element will on return point to the vector elements
        corresponding to the global Ids.
      \param wgts will on return point to a list of the weight or weights 
         associated with each element in the Ids list.  Weights are listed by 
         element by weight component.   NOT IMPLEMENTED YET
       \return The number of ids in the Ids list.
   */

  virtual size_t getMultiVectorView(int i, const gid_t *&Ids, 
     const scalar_t *&element, const scalar_t *&wgts) const = 0;

  /*! Given a new mapping of multi vertex elements to processes,
   *    create a new multi vertex with this mapping, and migrate
   *    the first vertex to it.  This is optional,
   *    This method only needs to be defined if you want to
   *    use it redistribute your multi vector. 
   *  TODO   documentation
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
