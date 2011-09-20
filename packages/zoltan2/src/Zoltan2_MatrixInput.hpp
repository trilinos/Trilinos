// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_MatrixInput.hpp

    \brief The abstract interface for a matrix input adapter.

    \author Siva Rajamanickam
*/


#ifndef _ZOLTAN2_MATRIXINPUT_HPP_
#define _ZOLTAN2_MATRIXINPUT_HPP_

#include <Zoltan2_InputAdapter.hpp>

namespace Zoltan2 {

/*! Zoltan2::MatrixInput
    \brief The MatrixInput is the abstract base class for matrix input adapters.

    The Matrix accessor methods defined here mimic those of Tpetra::CrsMatrix
    These public methods define the graph adapter interface to Zoltan2 models.
    TODO: It may be necessary to put a migration interface at this level.
*/

CONSISTENT_CLASS_TEMPLATE_LINE
  class MatrixInput : public InputAdapter{
private:

public:

  /*! Returns the global number of rows in the matrix.
   */
  virtual GNO getGlobalNumRows() const = 0;

  /*! Returns the global number columns in the matrix.
   */
  virtual GNO getGlobalNumCols() const = 0;

  /*! Returns the number rows on this process.
   */
  virtual LNO getNodeNumRows() const = 0;

  /*! Returns the number edges on this process.
   */
  virtual LNO getNodeNumCols() const = 0;

  //TODO: Add just the required functions.

};
  
  
}  //namespace Zoltan2
  
#endif
