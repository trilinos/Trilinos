#ifndef NumPyWrapper_h
#define NumPyWrapper_h

#ifndef numeric_include_h
#  include <numeric_include.h>
#endif

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
