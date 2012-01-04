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
  MatrixAdapterType,    /*!< matrix input */
  MeshAdapterType,    /*!< mesh input */
  GraphAdapterType,    /*!< graph input */
  CoordinateAdapterType,    /*!< coordinate input */
  VectorAdapterType,    /*!< vector input*/
  MultiVectorAdapterType,    /*!< multivector input*/
  IdentifierAdapterType,    /*!< plain identifier input*/
  XpetraCrsMatrixAdapterType  /*!< identify Xpetra adapters for better performance */
};


/*! \brief InputAdapter defines methods required by all InputAdapters
 */

class InputAdapter {
private:
public:

  /*! \brief Return type of adapter.
   */
  virtual enum InputAdapterType inputAdapterType()const = 0;

  /*! \brief Pure virtual destructor
   */
  virtual ~InputAdapter() {};

  /*! \brief Return a name that identifies the concrete adapter.
   */
  virtual std::string inputAdapterName() const = 0;

};
  
  
}  //namespace Zoltan2
  
#endif
