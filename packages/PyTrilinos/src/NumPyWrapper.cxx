#include "NumPyWrapper.h"

#include <iostream>

// Static initialization
NumPyWrapper         NumPyWrapper::m_singleton;

NumPyWrapper::NumPyWrapper()
{
  //std::cerr << "Calling import_array() to initialize Numeric C interface." << std::endl;
  import_array();
}

NumPyWrapper::~NumPyWrapper()
{
}

NumPyWrapper & NumPyWrapper::self()
{
  return m_singleton;
}

PyArrayObject * NumPyWrapper::pyContiguousArrayToDoubleArray(PyObject * p_pyObject                ,
                                                             double   * & numericArray            ,
                                                             int        & numericArrayNumDims     ,
                                                             int      * & p_numericArrayDimLengths,
                                                             int        & numericArrayTotalLength ,
                                                             int      * & p_numericArrayStrides    )
{
  PyArrayObject * array = contiguous_typed_array(p_pyObject, PyArray_DOUBLE, 0, 0);
  assert (0 != array && "Problem occurred.  Function returned a NULL pointer.");

  numericArrayNumDims      = (int)     array->nd        ;
  p_numericArrayDimLengths = (int *)   array->dimensions;
  numericArrayTotalLength  = PyArray_Size(p_pyObject)   ;
  numericArray             = (double *)array->data      ;
  p_numericArrayStrides    = (int *)   array->strides   ;

  return array;
}

PyArrayObject * NumPyWrapper::pyArrayToDoubleArray(PyObject * p_pyObject                ,
                                                   double   * & numericArray            ,
                                                   int        & numericArrayNumDims     ,
                                                   int      * & p_numericArrayDimLengths,
                                                   int        & numericArrayTotalLength ,
                                                   int      * & p_numericArrayStrides    )
{
  PyArrayObject * array = typed_array(p_pyObject, PyArray_DOUBLE, 0, 0);
  assert (0 != array && "Problem occurred.  Function returned a NULL pointer.");

  numericArrayNumDims      = (int)     array->nd        ;
  p_numericArrayDimLengths = (int *)   array->dimensions;
  numericArrayTotalLength  = PyArray_Size(p_pyObject)   ;
  numericArray             = (double *)array->data      ;
  p_numericArrayStrides    = (int *)   array->strides   ;

  return array;
}

void NumPyWrapper::LoadToExternalArray(const PyArrayObject      * p_pyArrayObject, 
                                       double                   * p_externalArray )
{
  const double * const p_pyArrayArray = (const double *)p_pyArrayObject->data;

  CopyToExternalArray(LOAD           ,
                      p_pyArrayObject, 
                      p_pyArrayArray , 
                      p_externalArray );

}
void NumPyWrapper::UnloadToExternalArray(PyArrayObject            * p_pyArrayObject, 
                                         const double             * p_externalArray )
{
  double       * const p_pyArrayArray = (double       *)p_pyArrayObject->data;
  
  CopyToExternalArray(UNLOAD         ,
                      p_pyArrayObject,
                      p_externalArray,
                      p_pyArrayArray  );
}

void NumPyWrapper::CopyToExternalArray(NumPyWrapper::CopyAction   copyAction     , 
                                       const PyArrayObject      * p_pyArrayObject,
                                       const double             * p_sourceArray  ,
                                       double                   * p_targetArray   )
{
  PyObject * const p_pyObject              = (PyObject *) p_pyArrayObject;
  assert (PyArray_Check(p_pyObject) && "Cast of PyArrayObject to PyObject failed.");

  const bool           isContiguous        = PyArray_ISCONTIGUOUS(p_pyArrayObject) ;
  const int            totalLength         = PyArray_Size(p_pyObject);
  
  if (isContiguous) {
    memcpy (p_targetArray, p_sourceArray , totalLength*sizeof(double));
  }
  else {
    const int            pyArrayNumDims      = (int)     p_pyArrayObject->nd         ;
    const int    * const p_pyArrayDimLengths = (int *)   p_pyArrayObject->dimensions ;
    const int    * const p_pyArrayStrides    = (int *)   p_pyArrayObject->strides    ;

    register       int   contiguousIndex     = 0;
    register       int   nonContiguousIndex  = 0;
    register       int   dimLength           = 0;
    register       int   result              = 0;
    register       int   remander            = 0;
    register       int   dimIndex            = 0;
    register       int   lastIndex           = 0;
    register const int   lastDim             = p_pyArrayDimLengths[pyArrayNumDims-1];
    register const int   lastStride          = p_pyArrayStrides[pyArrayNumDims-1]/sizeof(double);
    // Loop over all elements, striding by last dim's length
    for (contiguousIndex  = 0          ; 
         contiguousIndex  < totalLength; 
         contiguousIndex += lastDim     )
    {
      // Initialize non-contigous index calculation
      dimLength          = 0;
      remander           = 0;
      nonContiguousIndex = 0;
      result             = contiguousIndex;
      // Given contiguous index, calculate non-contiguous index
      switch (pyArrayNumDims)
      {
      case 1:
        dimLength           = p_pyArrayDimLengths      [0]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[0]/sizeof(double);
        break;
      case 2:
        dimLength           = p_pyArrayDimLengths      [1]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[1]/sizeof(double);
        result             /= dimLength                                  ;
        dimLength           = p_pyArrayDimLengths      [0]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[0]/sizeof(double);
        break;
      case 3:
        dimLength           = p_pyArrayDimLengths      [2]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[2]/sizeof(double);
        result             /= dimLength                                  ;
        dimLength           = p_pyArrayDimLengths      [1]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[1]/sizeof(double);
        result             /= dimLength                                  ;
        dimLength           = p_pyArrayDimLengths      [0]               ;
        remander            = result % dimLength                         ;
        nonContiguousIndex += remander*p_pyArrayStrides[0]/sizeof(double);
        break;
      default:
        for (dimIndex = 0; dimIndex < pyArrayNumDims; ++dimIndex)
        {
          dimLength  = p_pyArrayDimLengths[pyArrayNumDims-dimIndex-1];
          // remander is the array index for the
          // pyArrayNumDims-dimIndex-1 dimension
          remander   = result % dimLength;
          
          // Calculate nonContiguousIndex using NumPy array
          // stride-array-based formula
          nonContiguousIndex += remander*
            p_pyArrayStrides[pyArrayNumDims-dimIndex-1]/sizeof(double);
          result    /= dimLength;
        }
      }
      if (p_pyArrayStrides[pyArrayNumDims-1]/sizeof(double) == 1)
      {
        // Last stride is 1, so can use memcpy to fill
        if (copyAction == LOAD) {
          memcpy (&p_targetArray[contiguousIndex   ], 
                  &p_sourceArray[nonContiguousIndex], 
                  lastDim*sizeof(double)               );
        }
        else {
          memcpy (&p_targetArray[nonContiguousIndex], 
                  &p_sourceArray[contiguousIndex   ], 
                  lastDim*sizeof(double)               );
        }
      }
      else
      {
        const double * p_source    ;
        double       * p_target    ;
        register int   sourceStride;
        register int   targetStride;
        // Set ptr to begining of PyArray and external data using
        // non-contiguous and contiquous index.  
        // Set strides for final copy.
        if (copyAction == LOAD) {
          p_target     = &p_targetArray[contiguousIndex   ];
          p_source     = &p_sourceArray[nonContiguousIndex];
          targetStride = 1             ;
          sourceStride = lastStride    ;
        }
        else {
          p_target     = &p_targetArray[nonContiguousIndex];
          p_source     = &p_sourceArray[contiguousIndex   ];
          targetStride = lastStride    ;
          sourceStride = 1             ;
        }
        // Last stride not 1, so must use loop to fill
        // Loop over all elements in last dimension
        for (lastIndex = 0; lastIndex < lastDim; ++lastIndex)
        {
          // Copy data via assignment
          register const int targetIndex = lastIndex*targetStride;
          register const int sourceIndex = lastIndex*sourceStride;
          p_target[targetIndex]          = p_source[sourceIndex] ;
        } // Loop over last dimension
      } // Test if last dimension is contiguous
    } // Loop over contiguous index, striding by last dim's length
  } // Test if PyArray array is contiguous
}


//********************************************************************
// Tries to create a contiguous numeric array of type typecode from a 
// Python object. Works for list, tuples and numeric arrays.
//
// obj: Numeric array Python object
// typecode: data type PyArray_{ CHAR, UBYTE, SBYTE, SHORT, INT, LONG, FLOAT,
//                               DOUBLE, CFLOAT, CDOUBLE }
// expectnd: required number of dimensions. Used for checking. Ignored if <=0.
// expectdims: array of expected extends. Used for checking. Ignored if <=0.
//
// Raises ValueError exceptions if:
// - the PyArray_ContiguousFromObject fails
// - the array has a bad shape
// - the extent of a given dimension doesn't match the specified extent.
//********************************************************************
PyArrayObject * NumPyWrapper::contiguous_typed_array(PyObject *obj, int typecode,
                                                     int expectnd, int *expectdims)
{
  char buf[255];
  PyArrayObject *arr;

  

  // If the shape and type are OK, this function increments the reference
  // count and arr points to obj.
  if((arr = (PyArrayObject *)PyArray_ContiguousFromObject(obj,
                                                          typecode, 
                                                          0,
                                                          0)) == NULL)
  {
    sprintf(buf,"Failed to make a contiguous array of type %d\n", typecode);
    PyErr_SetString(PyExc_ValueError, buf);
    return NULL;
  }
  
  return array_check(arr, expectnd, expectdims, buf);
}

//********************************************************************
// Just like contiguous_typed_array() but gets a reference to the
// PyNum array.  This means you must use strides to access the elements
// of this possibly discontiguous array correctly.
// See contiguous_typed_array() for more information.
//********************************************************************
PyArrayObject * NumPyWrapper::typed_array(PyObject *obj, int   typecode  ,
                                          int expectnd , int * expectdims )
{
  char buf[255];
  PyArrayObject *arr;

  // If the shape and type are OK, this function increments the reference
  // count and arr points to obj.
  if((arr = (PyArrayObject *)PyArray_FromObject(obj,
                                                typecode, 
                                                0,
                                                0)) == NULL)
    {
      sprintf(buf,"Failed to get array of type %d\n", typecode);
      PyErr_SetString(PyExc_ValueError, buf);
      return NULL;
    }

  return array_check(arr, expectnd, expectdims, buf);
}


// Utility function for contiguous_typed_array() and typed_array().  
// See their documentation for more info.
PyArrayObject * NumPyWrapper::array_check(PyArrayObject * arr       ,
                                          int             expectnd  ,
                                          int           * expectdims,
                                          char          * buf        )
{
  // int i, numitems, itemsize;
  int i;

  if(expectnd>0)
    {
      if(arr->nd > expectnd + 1 || arr->nd < expectnd)
        {
          Py_DECREF((PyObject *)arr);
          PyErr_SetString(PyExc_ValueError,
                          "Array has wrong number of dimensions");
          return NULL;
        }
      if(arr->nd == expectnd + 1)
        {
          if(arr->dimensions[arr->nd - 1] != 1)
            {
              Py_DECREF((PyObject *)arr);
              PyErr_SetString(PyExc_ValueError,
                              "Array has wrong number of dimensions");
              return NULL;
            }
        }
      if(expectdims)
        {
          for(i = 0; i < expectnd; i++)
            if(expectdims[i]>0)
              if(expectdims[i] != arr->dimensions[i])
                {
                  Py_DECREF((PyObject *)arr);
                  sprintf(buf,"The extent of dimension %d is %d while %d was expected\n",
                          i, arr->dimensions[i], expectdims[i]);
                  PyErr_SetString(PyExc_ValueError, buf);
                  return NULL;
                }
                  
        }
    }

  return arr;
}


