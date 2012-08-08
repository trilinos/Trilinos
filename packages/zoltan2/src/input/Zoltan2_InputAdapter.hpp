// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_InputAdapter.hpp
    \brief Defines the InputAdapter interface.
*/

#ifndef _ZOLTAN2_INPUTADAPTER_HPP_
#define _ZOLTAN2_INPUTADAPTER_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_InputTraits.hpp>

namespace Zoltan2 {

/*! \brief An enum to identify general types of input adapters.
 *
 *  If you change this, update inputAdapterTypeName().
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

template <typename User>
  class InputAdapter {

private:

  typedef typename InputTraits<User>::scalar_t    scalar_t;

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
   *   Number of weights per object should be zero or greater.  If
   *   zero, then it is assumed that all objects are equally weighted.
   */ 
  virtual int getNumberOfWeightsPerObject() const = 0;

  /*! \brief Provide pointer to a weight array with stride.
   *    \param dim  the weight dimension, zero or greater
   *    \param wgt on return a pointer to the weights for this dimension
   *    \param stride on return, the value such that
   *       the \t nth weight should be found at <tt> wgt[n*stride] </tt>.
   *  \return the length of the \c wgt array, which should be at least
   *   equal to <tt> getLocalNumberOfObjects() * stride </tt>.
   */ 
  virtual size_t getObjectWeights(int dim, const scalar_t *&wgt, 
    int &stride) const = 0;

  /*! \brief Returns the name of the input adapter
   */
  static string inputAdapterTypeName(InputAdapterType iaType);
};

template <typename User>
  string InputAdapter<User>::inputAdapterTypeName(InputAdapterType iaType)
{
  string typeName;
  switch (iaType){
    case InvalidAdapterType:
      typeName = string("invalid");
      break;
    case IdentifierAdapterType:
      typeName = string("identifier");
      break;
    case VectorAdapterType:
      typeName = string("vector");
      break;
    case CoordinateAdapterType:
      typeName = string("coordinate");
      break;
    case GraphAdapterType:
      typeName = string("graph");
      break;
    case MeshAdapterType:
      typeName = string("mesh");
      break;
    case MatrixAdapterType:
      typeName = string("matrix");
      break;
    default:
      typeName = string("unknown");
      break;
  }

  return typeName;
}
  
  
}  //namespace Zoltan2
  
#endif
