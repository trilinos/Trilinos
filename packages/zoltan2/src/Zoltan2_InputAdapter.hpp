// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputAdapter.hpp

    \brief The abstract interface for an input adapter.
*/


#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

#include <Zoltan2_TemplateMacros.hpp>

namespace Zoltan2 {

/*! InputAdapter defines methods required by all InputAdapters
 */

enum InputAdapterType {
  InvalidAdapterType = 0,
  MatrixAdapterType,
  MeshAdapterType,
  GraphAdapterType,
  CoordAdapterType,
  IdAdapterType
};

class InputAdapter {
private:
public:
  virtual enum InputAdapterType adapterType() = 0;
protected:
};
  
  
}  //namespace Zoltan2
  
#endif
