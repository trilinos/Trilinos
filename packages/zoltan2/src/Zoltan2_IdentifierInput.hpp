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

#include <string>
#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::IdentifierInput
    \brief IdentifierInput defines the interface for input adapters 
           that represent a list of identifiers and weights.

    The template parameter is the the User's data structure.
    The user must have previously defined traits for their 
    data structure, supplying Zoltan2 mainly with their data types.
 
 TODO: weights
*/

template <typename User>
  class IdentifierInput : public InputAdapter {
private:

public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  enum InputAdapterType inputAdapterType() {return IdentifierAdapterType;}

  /*! Pure virtual destructor
   */
  virtual ~IdentifierInput() {};

  /*! Return the number of weights associated with each identifier.
   */
  virtual int getIdentifierWeightDim() const = 0;

  /*! Return the number of identifiers on this process.
   */
  virtual size_t getNumberOfIdentifiers() const = 0;

  /*! Sets pointers to this process' identifiers.
      \param Ids will on return point to the list of the global Ids for 
        this process.
      \param localIds can, optionally, on return point to a list of locally
        relevant values that the process will use to refer to the objects
        listed in the first list. If localIds is NULL and
        haveConsecutiveLocalIds is true, it is assumed that the
        global Ids are in local ID order.
      \param wgts will on return point to a list of the weight or weights 
         associated with each element in the Ids list.  Weights are listed by 
         element by weight component.   NOT IMPLEMENTED YET
       \return The number of ids in the Ids list.
   */

  virtual size_t getIdentifierView(const gid_t *&Ids,  const lid_t *&localIds,
     const scalar_t *&wgts) const = 0;

  /*! Given a new mapping of identifiers to partitions,
   *    migrate the identifiers to the new partitions.
   *    This method is optional, and only needs to be 
   *    defined if you want to redistribute your objects. 
   */
  size_t applyPartitioningSolution(User &in, User *&out,
    size_t numParts, size_t numIds,
    const gid_t *gid, const lid_t *lid, const size_t *partition)
  {
    return 0;
  } 
};
  
  
}  //namespace Zoltan2
  
#endif
