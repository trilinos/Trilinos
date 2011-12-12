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
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  enum InputAdapterType inputAdapterType() const {return IdentifierAdapterType;}

  /*! Pure virtual destructor
   */
  virtual ~IdentifierInput() {};

  /*! Return the number of weights associated with each identifier.
   */
  virtual int getNumWeights() const = 0;

  /*! Return the number of identifiers on this process.
   */
  virtual size_t getLocalNumIds() const = 0;

  /*! Sets pointers to this process' identifiers and optional weights.
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

  virtual size_t getIdList(const gid_t (*Ids), const scalar_t * (*weights),   
    int *strides) const = 0;


  /*! Given a new mapping of identifiers to partitions,
   *    migrate the identifiers to the new partitions.
   *    This method is optional, and only needs to be 
   *    defined if you want to redistribute your objects. 
   */
  size_t applyPartitioningSolution(User &in, User *&out,
    const PartitioningSolution<gid_t, lno_t> &solution)
  {
    return 0;
  } 
};
  
  
}  //namespace Zoltan2
  
#endif
