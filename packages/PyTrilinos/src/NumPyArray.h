#ifndef NumPyArray_h
#define NumPyArray_h

#ifndef numeric_include_h
#  include <numeric_include.h>
#endif

#include <iostream>

class NumPyArrayBase
{
public:
  NumPyArrayBase();
  virtual ~NumPyArrayBase () = 0;

  virtual std::ostream & print(std::ostream & stream) const = 0; 

  int                   getNumDims    () const;
  int                   getTotalLength() const;
  const int           * getDimLengths () const;
  double              * getDataArray  ()      ;
  const double        * getDataArray  () const;
  const int           * getStrides    () const;
  const PyArrayObject * getArrayObject() const;
  
  // Reports if wrapped array is contiguous
  bool isContiguous() const;

protected:
  // Print base class info.  To be called by implementation
  // of print.
  std::ostream & printBase(std::ostream & stream) const; 

  void setData(PyArrayObject * p_pyArrayObject       ,
               int             numPyArrayNumDims     ,
               int             numPyArrayTotalLength ,
               int           * p_numPyArrayDimLengths,
               int           * p_numPyArrayStrides   ,
               double        * p_dataArray            );

private:
  // Private and not implemented
  NumPyArrayBase(const NumPyArrayBase & a_ref);
  const NumPyArrayBase & operator = (const NumPyArrayBase & a_rhs);

private:
  PyArrayObject * mp_pyArrayObject       ;
  int             m_numPyArrayNumDims    ;
  int             m_numPyArrayTotalLength;
  int           * mp_numPyArrayDimLengths;
  int           * mp_numPyArrayStrides   ;
  double        * mp_dataArray           ;
};

class NumPyArrayContiguous: public NumPyArrayBase
{
public:
  NumPyArrayContiguous(PyObject * p_pyObject);
  virtual ~NumPyArrayContiguous ();

  virtual std::ostream & print(std::ostream & stream) const; 
  
private:
  // Private and not implemented
  NumPyArrayContiguous();
  NumPyArrayContiguous(const NumPyArrayContiguous & a_ref);
  const NumPyArrayContiguous & operator = (const NumPyArrayContiguous & a_rhs);
};

class NumPyArray: public NumPyArrayBase
{
public:
  NumPyArray(PyObject * p_pyObject);
  virtual ~NumPyArray ();

  virtual std::ostream & print(std::ostream & stream) const; 

private:
  // Private and not implemented
  NumPyArray();
  NumPyArray(const NumPyArray & a_ref);
  const NumPyArray & operator = (const NumPyArray & a_rhs);

};

#endif // NumPyArray_h
