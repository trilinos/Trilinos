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
  /*! Pure virtual destructor
   */
public:

  virtual ~Model() {};

  /*! An algorithm or third party library may
   *   require that global IDs be consecutive.
   *   It should be expected that this is not
   *   always as efficient as allowing arbitrary
   *   global IDs.
   *
   *   \param c  If true, a model will be created that has
   *      consecutive global ids.  If false, the most convenient
   *      global numbering will be used.
   */
  void setGlobalIdsMustBeConsecutive(bool c) {requireConsecutiveGlobalIds_=c;}

  /*! An algorithm or third party library may
   *   require that global IDs be consecutive.
   *
   *   \return  true, if a model will be created that has
   *      consecutive global ids.  false, if the most convenient
   *      global numbering will be used.
   */
  bool getGlobalIdsMustBeConsecutive()  { return requireConsecutiveGlobalIds_; }

private:

  bool requireConsecutiveGlobalIds_;
};

}

#endif
