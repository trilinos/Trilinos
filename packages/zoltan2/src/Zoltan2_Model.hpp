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
};

}

#endif
