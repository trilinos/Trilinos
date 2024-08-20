// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file NOX_LAPACK_Vector.H - Header file for NOX::LAPACK::Vector. */

#ifndef NOX_LAPACK_VECTOR_H
#define NOX_LAPACK_VECTOR_H

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H" // base class
#include <vector>
#include <ostream>

namespace NOX {

  //! NOX BLAS/LAPACK support
  namespace LAPACK {

  //! 1.0
  const double d_one = 1.0;
  //! -1.0
  const double d_mone = -1.0;
  //! 0.0
  const double d_zero = 0.0;
  //! 1
  const int i_one = 1;
  //! 0
  const int i_zero = 0;

    //! Implementation of NOX::Abstract::Vector for STL std::vector<double> (using LAPACK for some computations)
    class Vector : public NOX::Abstract::Vector {

    public:

      //! Construct an empty vector
      Vector();

      //! Construct a zero vector of length n
      Vector(int n);

      //! Construct a vector of length n from given array
      Vector(int n, double *v);

      //! Copy constructor
      Vector(const NOX::LAPACK::Vector& source,
         NOX::CopyType type = NOX::DeepCopy);

      //! Destruct Vector.
      ~Vector();

      //@{ \name Initialization methods.

      NOX::Abstract::Vector& init(double value);

      //! Initialize every element of this vector with random values
      virtual NOX::Abstract::Vector& random(bool useSeed = false, int seed =  1);

      //! Replace this vector with STL double vector y
      NOX::Abstract::Vector& operator=(const std::vector<double>& y);

      // derived
      NOX::Abstract::Vector& operator=(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);

      // derived
      NOX::Abstract::Vector& abs(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

      // derived
      NOX::Abstract::Vector& reciprocal(const NOX::LAPACK::Vector& y);
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

      //@}

      //@{ \name Update methods.

      // derived
      NOX::Abstract::Vector& scale(double gamma);

      // derived
      NOX::Abstract::Vector& scale(const NOX::LAPACK::Vector& a);
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::LAPACK::Vector& a, double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a, double gamma = 0.0);

      // derived
      NOX::Abstract::Vector& update(double alpha, const NOX::LAPACK::Vector& a,
                   double beta, const NOX::LAPACK::Vector& b,
                   double gamma = 0.0);
      NOX::Abstract::Vector& update(double alpha, const NOX::Abstract::Vector& a,
                   double beta, const NOX::Abstract::Vector& b,
                   double gamma = 0.0);

      //@}

      //@{ \name Creating new Vectors.

      // derived
      Teuchos::RCP<NOX::Abstract::Vector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //@}

      //@{ \name Norms.

      // derived
      double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;

      // derived
      double norm(const NOX::LAPACK::Vector& weights) const;
      double norm(const NOX::Abstract::Vector& weights) const;

      //@}

      //@{ \name Inner products

      // derived
      double innerProduct(const NOX::LAPACK::Vector& y) const;
      double innerProduct(const NOX::Abstract::Vector& y) const;

      //@}

      // derived
      NOX::size_type length() const;

      //! Return the i-th element
      double& operator() (int i);

      //! Return the i-th element (const version)
      const double& operator() (int i) const;

      //! Prints out the vector to the specified stream.
      /*!
    For example, a vector would appear as
    \f[ \left[ \; 0.1 \; 2.34 \; 5 \; \right] \f]
    It will be all on one line, with a single space between each entry, bracketed on either side.
      */
      std::ostream& leftshift(std::ostream& stream) const;

      // derived
      void print(std::ostream& stream) const;

    private:

      //! Return the i-th element
      double& operator[] (int i);

      //! Return the i-th element
      const double& operator[] (int i) const;

      //! The used length of vector
      int n;

      //! The vector owned by this object
      std::vector<double> x;
    };

    //! Function for printing
    std::ostream& operator<<(std::ostream& stream, const NOX::LAPACK::Vector& v);

  } // namespace LAPACK
} // namespace NOX




#endif
