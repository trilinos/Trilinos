// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IdentifierAdapter.hpp

    \brief The abstract interface for an simple identifier input.
*/


#ifndef _ZOLTAN2_IDENTIFIERADAPTER_HPP_
#define _ZOLTAN2_IDENTIFIERADAPTER_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierTraits.hpp>

namespace Zoltan2 {

/*! \brief IdentifierAdapter defines methods required by all input adapters
 *    for simple identifer input.
 *
 *  This interface is sufficient when the objects to be partitioned
 *  or ordered are simple global identifiers with weights.
 */

template <typename User>
class IdentifierAdapter : InputAdapter {
private:
public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;

  /*! \brief Return type of adapter.
   */
  enum InputAdapterType inputAdapterType() { return IdAdapterType};

  /*! \brief Pure virtual destructor
   */
  virtual ~IdentifierAdapter() {};

  /*! \brief Return the number of identifiers on the local process.
   */
  virtual size_t getLocalNumIds() const = 0;

  /*! \brief Return the local identifiers.
   *   \param gids on return should point to the global identifiers
   *   \param if the caller is using local identifiers to reference
   *      the global identifiers, and the local identifiers are not
   *      consecutive integers, then a pointer to the local identifiers 
   *      corresponding to the global identifiers must be set here.
       \return the number of identifiers in the lists.
   */
  virtual size_t getIdList(const gid_t *&gids, const lid_t *&lids) const=0;

  /*! \brief Return the number of weights per identifier.
   *
   *   If the number of weights is zero, the Zoltan2 algorithm will
   *   assume that all identifiers should be weighted equally.
   */
  virtual int getNumWeights();

  /*! \brief Return one set of weights.
   *   \ param n  The index for the weights to return, ranging from
   *              zero to getNumWeights() - 1.
   *   \ param weights  on return should point to the weights for
   *     index n corresponding to the Identifiers returned in getIdList. 
   */
  virtual void getWeights(int n, scalar_t *&weights);
};
  
  
}  //namespace Zoltan2
  
#endif
