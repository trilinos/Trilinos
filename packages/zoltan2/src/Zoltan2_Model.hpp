// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Model.hpp

    \brief The abstract interface for a computational model.
*/


#ifndef _ZOLTAN2_MODEL_HPP_
#define _ZOLTAN2_MODEL_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_IdentifierMap.hpp>

namespace Zoltan2 {

/*! Zoltan2::Model
    to be called asynchronously.
*/
enum ModelType {
  InvalidModel = 0,
  HypergraphModelType,
  GraphModelType,
  GeometryModelType,
  IdModelType
};


template <typename Adapter>
class Model
{
public:

   /*! Zoltan2 identifier types, user identifier types
    */
  typedef typename Adapter::lno_t    lno_t;
  typedef typename Adapter::gno_t    gno_t;
  typedef typename Adapter::lid_t    lid_t;
  typedef typename Adapter::gid_t    gid_t;
  typedef IdentifierMap<lid_t, gid_t, lno_t, gno_t> idmap_t;

  /*! Pure virtual destructor
   */
  virtual ~Model() {};

  /*! Constructor
   */
  Model() : requireConsecutiveIds_(false), idMap_() {}

  /*! An algorithm or TPL may require consecutive global IDs.
   *
   * \param c  If true, a model will be created that has
   *   consecutive global ids.  If false, the most convenient
   *   global numbering will be used.
   */

  void setRequireConsecutiveIds(bool c) {requireConsecutiveIds_=c;}

  /*! An algorithm or TPL may require consecutive global IDs.
   *
   *   \return  true, if a model will be created that has
   *      consecutive global ids.  false, if the most convenient
   *      global numbering will be used.
   */
  bool getRequireConsecutiveIds()  { return requireConsecutiveIds_; }

  /*! Return the IdentifierMap used by the model.
   *
   *  Every model must have an IdentifierMap, whether it needs for mapping 
   *  or not. The Map can simply indicate that Zoltan2 global numbers are identical
   *  to the application's global IDs.
   */
  template <typename lid_t, typename gid_t, typename lno_t, typename gno_t>
    const RCP<const Zoltan2::IdentifierMap<lid_t, gid_t, lno_t, gno_t> > &getIdentifierMap()
    {
      return idMap_;
    }

protected:

  /*! Set the IdentifierMap used by the model.
   */
  template <typename lid_t, typename gid_t, typename lno_t, typename gno_t>
    void setIdentifierMap(RCP<const idmap_t> &map)
    {
      idMap_ = map; 
    }

private:

  bool requireConsecutiveIds_;
  RCP<const idmap_t> idMap_;
};

}

#endif
