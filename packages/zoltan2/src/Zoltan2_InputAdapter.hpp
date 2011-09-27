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

#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! InputAdapter defines methods required by all InputAdapters
 *
 *  About local IDs:
 *    Applications are required to supply unique global IDs for
 *    their objects, such as vertices, coordinates, matrix rows,
 *    and so on.
 *
 *    Local IDs are optional.  Local IDs are symbols which are
 *    meaningful to the application and that reference the object.
 *    A local ID may be an array index or a pointer for example.
 *
 *    The impact on input adapter, models and solutions is this:
 *       1. If local IDs are supplied, they must appear in the
 *          solution if the application requests a solution.
 *       2. The "set" methods in the input adapter must accept
 *          local IDs and global IDs, and the application may
 *          use either or both.
 *       3. If local IDs are available, the model should query
 *          the input adapter using local IDs, since normally this 
 *          is more efficient.
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
  // Return enumerated InputAdapterType for the input adapter.
  // This function is implemented in the MatrixAdapter, GraphAdapter,
  // MeshAdapter, CoordAdapter and IdAdapter subclasses.
  // Users do not have to implement this function for their adapters
  // as long as they inherit from one of the subclasses (which they must).
  virtual enum InputAdapterType inputAdapterType() = 0;

  /*! Return a name that identifies the concrete adapter.
   *  Useful for debugging.
   */
  virtual std::string inputAdapterName() const = 0;

protected:
};
  
  
}  //namespace Zoltan2
  
#endif
