// @HEADER
// ***********************************************************************
//
//                 PyTrilinos: Rapid Prototyping Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef NumPyWrapper_h
#define NumPyWrapper_h

#include "numpy_include.h"

// Singleton class that wraps the various Numeric C interface.
// Important step is the calling of the function import_array() in the constructor.
// This class is a singleton that gets instantiated before main() is called (because
// this class has a static attribute of itself) so that import_array() will be
// called before anything else happens.

class NumPyWrapper
{
public:
  enum CopyAction
    {
      LOAD   = 0,
      UNLOAD
    };
  
  static PyArrayObject * contiguous_typed_array(PyObject *obj, int   typecode  ,
                                                int expectnd , int * expectdims );

  static PyArrayObject * typed_array(PyObject *obj, int   typecode  ,
                                     int expectnd , int * expectdims );
  
  static PyArrayObject * pyContiguousArrayToDoubleArray(PyObject * p_pyObject                ,
                                                        double   * & numericArray            ,
                                                        int        & numericArrayNumDims     ,
                                                        int      * & p_numericArrayDimLengths,
                                                        int        & numericArrayTotalLength ,
                                                        int      * & p_numericArrayStrides   );

  static PyArrayObject * pyArrayToDoubleArray(PyObject * p_pyObject                ,
                                              double   * & numericArray            ,
                                              int        & numericArrayNumDims     ,
                                              int      * & p_numericArrayDimLengths,
                                              int        & numericArrayTotalLength ,
                                              int      * & p_numericArrayStrides    );

  static NumPyWrapper & self();


  static void LoadToExternalArray  (const PyArrayObject      * p_pyArrayObject, 
                                    double                   * p_externalArray );

  static void UnloadToExternalArray(PyArrayObject            * p_pyArrayObject, 
                                    const double             * p_externalArray );

  static void CopyToExternalArray  (NumPyWrapper::CopyAction   copyAction     , 
                                    const PyArrayObject      * p_pyArrayObject,
                                    const double             * p_sourceArray  ,
                                    double                   * p_targetArray   );

private:
  static PyArrayObject * array_check(  PyArrayObject * arr       ,
                                       int             expectnd  ,
                                       int           * expectdims,
                                       char          * buf        );

protected:
  // These are protected instead of private to keep
  // compilers happy.
  ~NumPyWrapper();
  NumPyWrapper ();

private:
  NumPyWrapper(const NumPyWrapper & a_ref);
  const NumPyWrapper & operator = (const NumPyWrapper & a_rhs);

private:
  // The singleton, i.e. the only instance of this object is this
  // attribute
  static NumPyWrapper m_singleton;
  
};

#endif // NumPyWrapper_h
