#include "NumPyArray.h"
#include "NumPyWrapper.h"

#include <iostream>
#include <string>
#include <sstream>

// Utlity Function
template<typename T>
std::string ToStringFrom(const T& t) {
  std::ostringstream s;
  s << t;
  return s.str();
}


// NumPyArrayBase -------------------------------------------------------------
NumPyArrayBase::NumPyArrayBase() :
  mp_pyArrayObject       (0),
  m_numPyArrayNumDims    (0),
  m_numPyArrayTotalLength(0),
  mp_numPyArrayDimLengths(0),
  mp_numPyArrayStrides   (0),
  mp_dataArray           (0)
{
  // std::cout << "In NumPyArrayBase::NumPyArrayBase()" << std::endl;
}

void NumPyArrayBase::setData(PyArrayObject * p_pyArrayObject       ,
                             int             numPyArrayNumDims     ,
                             int             numPyArrayTotalLength ,
                             int           * p_numPyArrayDimLengths,
                             int           * p_numPyArrayStrides   ,
                             double        * p_dataArray            )
{
  assert (0 != p_dataArray            && "Null pointer sent to function.");
  assert (0 != p_numPyArrayDimLengths && "Null pointer sent to function.");
  assert (0 != p_numPyArrayStrides    && "Null pointer sent to function.");
  assert (0 != p_pyArrayObject        && "Null pointer sent to function.");

  mp_pyArrayObject        = p_pyArrayObject       ;
  m_numPyArrayNumDims     = numPyArrayNumDims     ;
  m_numPyArrayTotalLength = numPyArrayTotalLength ;
  mp_numPyArrayDimLengths = p_numPyArrayDimLengths;
  mp_numPyArrayStrides    = p_numPyArrayStrides   ;
  mp_dataArray            = p_dataArray           ;
}

NumPyArrayBase::~NumPyArrayBase ()
{
  Py_DECREF((PyObject *) mp_pyArrayObject);
  
  mp_pyArrayObject        = 0;
  m_numPyArrayNumDims     = 0;
  m_numPyArrayTotalLength = 0;
  mp_numPyArrayDimLengths = 0;
  mp_numPyArrayStrides    = 0;
  mp_dataArray            = 0;
}

int NumPyArrayBase::getNumDims () const
{
  return m_numPyArrayNumDims;
}

int NumPyArrayBase::getTotalLength() const
{
  return m_numPyArrayTotalLength;
}

const int * NumPyArrayBase::getDimLengths () const
{
  return mp_numPyArrayDimLengths;
}

double * NumPyArrayBase::getDataArray ()
{
  return mp_dataArray;
}

const double * NumPyArrayBase::getDataArray () const
{
  return mp_dataArray;
}

const PyArrayObject * NumPyArrayBase::getArrayObject() const
{
  return mp_pyArrayObject;
}

PyArrayObject * NumPyArrayBase::getArrayObject()
{
  return mp_pyArrayObject;
}

bool NumPyArrayBase::isContiguous() const
{
  assert (0 != mp_pyArrayObject &&
          "No PyArrayObject object present in this object"
          "This should be impossible"                     );

  return PyArray_ISCONTIGUOUS(mp_pyArrayObject);
}

std::ostream & NumPyArrayBase::printBase(std::ostream & stream) const
{
  assert (0 != mp_pyArrayObject && 
          "Object in uninitialized state.  This should be impossible");

  register int   index  ;
  char         * p_comma;

  // Print array info
  stream << "Total length               = " << m_numPyArrayTotalLength << std::endl
         << "Number of dimensions       = " << m_numPyArrayNumDims     << std::endl
         << "Length in each dimension   = (";
  p_comma = "";
  for (index = 0; index < m_numPyArrayNumDims; ++index) {
    stream << p_comma << mp_numPyArrayDimLengths[index];
    p_comma = ", ";
  }
  stream << ")" << std::endl
         << "Is contiguous              = " << (isContiguous() ? "True" : "False") << std::endl
         << "Strides for each dimension = (";
  p_comma = "";
  for (index = 0; index < m_numPyArrayNumDims; ++index) {
    stream << p_comma << mp_numPyArrayStrides[index]/sizeof(double);
    p_comma = ", ";
  }
  stream << ")" << std::endl;

  // Print array data

  register int i, j     ;
  std::string  tmpString;
  for (i = 0; i < m_numPyArrayTotalLength; ++i) {
    int length   = 0;
    int result   = i;
    int remander = 0;
    index        = 0;
    std::string indexString = "";
    p_comma      = "";
    for (j = 0; j < m_numPyArrayNumDims; ++j) {
      length       = mp_numPyArrayDimLengths[m_numPyArrayNumDims-j-1];
      remander     = result   % length;
      tmpString    = ToStringFrom(remander);
      tmpString   += p_comma;
      indexString  = tmpString + indexString;
      p_comma      = ",";
      index += remander*mp_numPyArrayStrides[m_numPyArrayNumDims-j-1]/sizeof(double);
      result   /= length;
    }
    stream << "[" << indexString << "] (" << i << ") = " 
           << mp_dataArray[index] << std::endl;
  }
  
  return stream;
}

const int * NumPyArrayBase::getStrides () const
{
  return mp_numPyArrayStrides;
}

// NumPyArray -----------------------------------------------------------------

NumPyArray::NumPyArray(PyObject * p_pyObject) :
  NumPyArrayBase()  
{
  assert (0 != p_pyObject && "NULL pointer passed to constructor");

  int             numPyArrayNumDims         ;
  int             numPyArrayTotalLength     ; 
  int           * p_numPyArrayDimLengths = 0;
  int           * p_numPyArrayStrides    = 0;
  double        * p_dataArray            = 0;

  PyArrayObject * p_pyArrayObject = NumPyWrapper::
    pyArrayToDoubleArray(p_pyObject            ,
                         p_dataArray           ,
                         numPyArrayNumDims     ,
                         p_numPyArrayDimLengths,
                         numPyArrayTotalLength ,
                         p_numPyArrayStrides    );
  
  assert (0 != p_pyArrayObject        && "Function should have set this pointer");
  assert (0 != p_dataArray            && "Function should have set this pointer");
  assert (0 != p_numPyArrayDimLengths && "Function should have set this pointer");
  assert (0 != p_numPyArrayStrides    && "Function should have set this pointer");
  
  setData(p_pyArrayObject       ,
          numPyArrayNumDims     ,
          numPyArrayTotalLength ,
          p_numPyArrayDimLengths,
          p_numPyArrayStrides   ,
          p_dataArray            );
}

NumPyArray::~NumPyArray ()
{
  // Does nothing here yet
}

std::ostream & NumPyArray::print(std::ostream & stream) const
{
  return printBase(stream);
}

// NumPyArrayContiguous -------------------------------------------------------

NumPyArrayContiguous::NumPyArrayContiguous(PyObject * p_pyObject) :
  NumPyArrayBase()
{
  // std::cout << "In NumPyArrayContiguous::NumPyArrayContiguous(...)" << std::endl;
  assert (0 != p_pyObject && "NULL pointer passed to constructor");

  int             numPyArrayNumDims         ;
  int             numPyArrayTotalLength     ; 
  int           * p_numPyArrayDimLengths = 0;
  int           * p_numPyArrayStrides    = 0;
  double        * p_dataArray            = 0;

  PyArrayObject * p_pyArrayObject = NumPyWrapper::
    pyContiguousArrayToDoubleArray(p_pyObject            ,
                                   p_dataArray           ,
                                   numPyArrayNumDims     ,
                                   p_numPyArrayDimLengths, 
                                   numPyArrayTotalLength ,
                                   p_numPyArrayStrides    );
  
  assert (0 != p_numPyArrayDimLengths && "Function should have set this pointer");
  assert (0 != p_dataArray            && "Function should have set this pointer");
  assert (0 != p_pyArrayObject        && "Function should have set this pointer");
  assert (0 != p_numPyArrayStrides    && "Function should have set this pointer");
  
  setData(p_pyArrayObject       ,
          numPyArrayNumDims     ,
          numPyArrayTotalLength ,
          p_numPyArrayDimLengths,
          p_numPyArrayStrides   ,
          p_dataArray            );
}

NumPyArrayContiguous::~NumPyArrayContiguous ()
{
  // Does nothing here yet
}

std::ostream & NumPyArrayContiguous::print(std::ostream & stream) const
{
  return printBase(stream);
}
