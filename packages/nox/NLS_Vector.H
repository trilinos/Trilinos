#ifndef _NLS_VECTOR_H
#define _NLS_VECTOR_H
 
//! NLSPACK pure abstract vector interface for vectors that are used
//! by the NLS solver. We require only limited functionality for the
//! vectors.

class NLS_Vector {

protected:

  //! NLS_Vector constructor
  //NLS_Vector();

public:

  //! NLS_Vector deconstructor
  ~NLS_Vector(){};

  //@{ \name Value initialization methods
  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
  virtual void random() = 0;

  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

  */
  virtual void putScalar (double scalar) = 0;

  //@{ \name Vector-by-vector methods
  //! Update multi-vector values with scaled values of A, \e this = scalar*\e this + scalarA*A.
  /*!
    \param In
           scalarA - scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           scalar - scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

  */
  virtual void update(double scalarA, const NLS_Vector& A, double scalar) = 0;

  //! Scale vector by alpha
  /*! Let x represent this vector. 
    Compute x = alpha * x. 
    Return x. */
  virtual void scale(double alpha) = 0;

  //! Dot product with y
  /*!  
    Compute value = <\e this,y> = sum_i (\e this[i] * y[i]).  
    Return value */
  virtual double dot(const NLS_Vector& y) const = 0;
  
  //! Infinity-norm
  /*! Let x represent this vector. 
    Compute value = x = max(|x|).
    Return value.
  */
  virtual double normInf() const = 0;

  //! 1-norm
  /*! Let x represent this vector. 
    Compute value = 
    Return value.
  */
  virtual double norm1() const = 0;

  //! 2-norm
  /*! Let x represent this vector. 
    Compute value = ||x||_2 = sqrt(sum_i x[i]^2).
    Return value.
  */
  virtual double norm2() const = 0;

  //! Weighted 2-norm
  /*! Let x represent this vector. 
    Compute value = ||x||_W = sqrt(sum_i (weights[i] * (x[i]^2)))??
    Return value.
  */
  virtual double norm2(const NLS_Vector& weights) const = 0;

  //! Puts element-wise absolute values of input Multi-vector in target.
  /*!
    \param In
           y - Input Vector.
    \param Out
           \e this will contain the absolute values of the entries of y.

    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for y and \e this.
  */
  virtual void Abs(const NLS_Vector& y) = 0;

  //! Compute minimum value of each vector in multi-vector.
  /*!
    \return minimum value of vector.
  */
  virtual double minValue() const = 0;

  //! Compute maximum value of vector.
  /*!
    \return maximum value of vector.
  */
  virtual double maxValue() const = 0;

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \return mean (average) value of vector.
  */
  virtual double meanValue() const = 0;

};

#endif
