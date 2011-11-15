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
#include <Zoltan2_InputTraits.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of input adapters.
 */
enum InputAdapterType {
  InvalidAdapterType = 0,    /*!< unused value */
  MatrixAdapterType,    /*!< matrix input */
  MeshAdapterType,    /*!< mesh input */
  GraphAdapterType,    /*!< graph input */
  CoordAdapterType,    /*!< coordinate input */
  VectorAdapterType,    /*!< vector input*/
  MultiVectorAdapterType,    /*!< multivector input*/
  IdAdapterType,    /*!< plain identifier input*/
  XpetraCrsMatrixAdapterType  /*!< identify Xpetra adapters for better performance */
};


/*! \brief InputAdapter defines methods required by all InputAdapters
 *
 *  Regarding local IDs:
 *  - Applications are required to supply unique global IDs for
 *    their objects, which may be vertices, coordinates, matrix rows,
 *    and so on.
 *
 *  - Local IDs are optional.  Local IDs are symbols which are
 *    meaningful to the application and which reference the object.
 *    A local ID may be an array index or a pointer, for example.
 *
 *  - The impact on input adapter, models and solutions is this:
 *       - If local IDs are supplied, they must appear in the
 *          solution if the application requests a solution.
 *       - The "set" methods in the input adapter must accept
 *          local IDs and global IDs, and the application may
 *          use either or both.
 *       - If local IDs are available and are consecutive integers, 
            the model should query
 *          the input adapter using local IDs, since normally this 
 *          is more efficient.
 */

class InputAdapter {
private:
public:

  /*! \brief Return type of adapter.
   */
  virtual enum InputAdapterType inputAdapterType() = 0;

  /*! \brief Pure virtual destructor
   */
  virtual ~InputAdapter() {};

  /*! \brief Return a name that identifies the concrete adapter.
   */
  virtual std::string inputAdapterName() const = 0;

  /*! \brief Returns true if input adapter uses local Ids for objects.  */
  virtual bool haveLocalIds() const = 0;

  /*! \brief Return true if local Ids are consecutive integral
   *   values and supply the base.  
   *
   *  Providing this information
   *   can save memory, making local ID lists unneccesary.
   */
  virtual bool haveConsecutiveLocalIds(size_t &base) const = 0;

};
  
  
}  //namespace Zoltan2
  
#endif
