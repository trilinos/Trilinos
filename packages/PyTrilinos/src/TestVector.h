#ifndef TestVector_h
#define TestVector_h

#include <vector>

#include <Python.h>

class TestVector
{
public:

  TestVector (int length);
  ~TestVector();

  const std::vector<double> & getVector() const {
    return m_vector;
  }

  void setVector(int i, double val) {
    m_vector[i] = val;
  }

  void fromNumericArray(PyObject * p_pyObject);
  void testit(PyObject * p_pyObject);
  void write() const;

#ifdef SWIG
  // Some helper functions
  %extend {
    double get (int i) {
      return self->getVector()[i];
    }
    void set(int i, double val) {
      self->setVector(i,val);
    }
  }

  //%include <NumericTypeMaps.i>
#endif

private:
  // Private and not implemented so never can be called
  TestVector ();
  TestVector(const TestVector & a_ref);
  const TestVector & operator = (const TestVector & a_rhs);

  std::vector<double> m_vector;
};

#endif // TestVector_h
