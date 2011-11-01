// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_VectorInput.hpp

    \brief The abstract interface for an input adapter representing a 
    distributed vector. 
*/


#ifndef _ZOLTAN2_VECTORINPUT_HPP_
#define _ZOLTAN2_VECTORINPUT_HPP_

#include <string>
#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::VectorInput
    \brief VectorInput defines the interface for input adapters for vectors.

    lid_t: the type for the application's local Ids
    gid_t: the type for the application's global Ids
    lno_t: the integral type that Zoltan2 will use for local counters.
    gno_t: the integral type that Zoltan2 will use for the global 
      counts and identifiers.  It needs to be large enough for the
      problem's number of objects.

    The template parameter is the the User's vector data structure.
    The user must have previously defined traits for their vector 
    data structure, supplying Zoltan2 mainly with their data types.
 
 TODO: weights and coordinates
*/

template <typename User>
  class VectorInput : public InputAdapter<User> {
private:

public:

  typedef typename InputAdapter<User>::scalar_t scalar_t;
  typedef typename InputAdapter<User>::lno_t    lno_t;
  typedef typename InputAdapter<User>::gno_t    gno_t;
  typedef typename InputAdapter<User>::lid_t    lid_t;
  typedef typename InputAdapter<User>::gid_t    gid_t;
  typedef typename InputAdapter<User>::node_t   node_t;

  enum InputAdapterType inputAdapterType() {return VectorAdapterType;}

  /*! Pure virtual destructor
   */
  virtual ~VectorInput() {};

  /*! Return the length of the portion of the vector on this process.
   */
  virtual size_t getLocalLength() const = 0;

  /*! Return the global length of the vector on this process.
   */
  virtual size_t getGlobalLength() const = 0;

  /** TODO getStride - what exactly does this mean
   */

  /*! Sets pointers to this process' vertex elements.
      \param Ids will on return point to the list of the global Ids for 
        each element on this process.
      \param localIds can, optionally, on return point to a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list. If localIds is NULL and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local ID order.
      \param element will on return point to the vector elements
        corresponding to the global Ids.
      \param wgts will on return point to a list of the weight or weights 
         associated with each element in the Ids list.  Weights are listed by 
         element by weight component.   NOT IMPLEMENTED YET
       \return The number of ids in the Ids list.
   */

  virtual lno_t getVectorView(gid_t *&Ids,  *&localIds,
     scalar_t *&element, scalar_t *&wgts) const = 0;

  /*! Given a new mapping of vertex elements to processes,
   *    create a new vertex with this mapping, and migrate
   *    the first vertex to it.  This is optional,
   *    This method only needs to be defined if you want to
   *    use it redistribute your vector. 
   *  TODO   documentation
   */
  int applyPartitioningSolution(User &in, User *&out,
    int numIds, int numParts, gid_t *gid,  *lid, int *partition) 
  {
    return 0;
  } 
};
  
  
}  //namespace Zoltan2
  
#endif
