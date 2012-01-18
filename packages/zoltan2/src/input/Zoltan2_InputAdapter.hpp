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

/*! \brief An enum to identify general types of input adapters.
 */
enum InputAdapterType {
  InvalidAdapterType = 0,    /*!< unused value */
  IdentifierAdapterType,    /*!< plain identifier input, just a list of Ids*/
  VectorAdapterType,    /*!< vector input*/
  CoordinateAdapterType,    /*!< coordinate input */
  GraphAdapterType,    /*!< graph input */
  MeshAdapterType,    /*!< mesh input */
  MatrixAdapterType    /*!< matrix input */
};


/*! \brief InputAdapter defines methods required by all InputAdapters

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.
 */

class InputAdapter {
private:
public:

  /*! \brief Returns the type of adapter.
   */
  virtual enum InputAdapterType inputAdapterType()const = 0;

  /*! \brief Pure virtual destructor
   */
  virtual ~InputAdapter() {};

  /*! \brief Returns a descriptive name that identifies the concrete adapter.
   */
  virtual std::string inputAdapterName() const = 0;

};
  
  
}  //namespace Zoltan2
  
#endif
