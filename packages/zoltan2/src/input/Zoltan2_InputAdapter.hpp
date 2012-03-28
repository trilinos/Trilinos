// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_InputAdapter.hpp
    \brief Defines the InputAdapter interface.
*/

#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

#include <Zoltan2_Standards.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of input adapters.
 */
enum InputAdapterType {
  InvalidAdapterType = 0,    /*!< \brief unused value */
  IdentifierAdapterType,    /*!< \brief plain identifier input, just a list of Ids*/
  VectorAdapterType,    /*!< \brief vector input*/
  CoordinateAdapterType,    /*!< \brief coordinate input */
  GraphAdapterType,    /*!< \brief graph input */
  MeshAdapterType,    /*!< \brief mesh input */
  MatrixAdapterType    /*!< \brief matrix input */
};


/*! \brief InputAdapter defines methods required by all InputAdapters

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    \todo Add add a MeshInput adapter
 */

class InputAdapter {
private:
public:

  /*! \brief Returns the type of adapter.
   */
  virtual enum InputAdapterType inputAdapterType()const = 0;

  /*! \brief Desstructor
   */
  virtual ~InputAdapter() {};

  /*! \brief Returns a descriptive name that identifies the concrete adapter.
   */
  virtual string inputAdapterName() const = 0;

  /*! \brief Returns the number of objects in the input.
   *
   *  Objects may be coordinates, graph vertices, matrix rows, etc.
   *  They are the objects to be partitioned, ordered, or colored.
   */
  virtual size_t getLocalNumberOfObjects() const = 0;

  /*! \brief Returns the number of weights per object.
   */ 
  virtual int getNumberOfWeightsPerObject() const = 0;
};
  
  
}  //namespace Zoltan2
  
#endif
