// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

/*! \file NOX_Thyra_Vector.H - Header file for NOX::Thyra::Vector. */

#ifndef NOX_THYRA_VECTOR_H
#define NOX_THYRA_VECTOR_H

#include "NOX_Abstract_Vector.H" // base class
#include "NOX_Abstract_ImplicitWeighting.H" // base class
#include "Thyra_VectorBase.hpp"
#include "Teuchos_RCP.hpp"

#ifndef TEUCHOS_DISABLE_ALL_TIMERS
// Define this to see selected timers
#define ENABLE_NOX_THYRA_TIMERS
#endif

#ifdef ENABLE_NOX_THYRA_TIMERS
#include "Teuchos_TimeMonitor.hpp"
#endif

namespace NOX {

  //! NOX Thyra support
  namespace Thyra {

    //! Implementation of NOX::Thyra::Vector
    class Vector : public NOX::Abstract::Vector,
           public NOX::Abstract::ImplicitWeighting {

    public:

      //! Construct a vector as a view of a given Thyra vector
      /*! Creates a view and does not allocate a new vector */
      Vector(const Teuchos::RCP< ::Thyra::VectorBase<double> >& src);

      //! Construct a vector from a given Thyra vector
      /*! Allocates an entirely new vector */
      Vector(const ::Thyra::VectorBase<double>& source);

      //! Copy constructor
      Vector(const NOX::Thyra::Vector& source,
         NOX::CopyType type = NOX::DeepCopy);

      //! Destruct Vector.
      ~Vector();

     //@{ \name Access to underlying Thyra vector.

      //! Get reference to underlying Thyra vector.
      virtual ::Thyra::VectorBase<double>& getThyraVector();

      //! Get const reference to underlying Thyra vector.
      virtual const ::Thyra::VectorBase<double>& getThyraVector() const;

      //! Get RCP to underlying Thyra vector.
      virtual Teuchos::RCP< ::Thyra::VectorBase<double> > getThyraRCPVector();

      //! Get RCP to underlying Thyra vector.
      virtual Teuchos::RCP<const ::Thyra::VectorBase<double> > getThyraRCPVector() const;

      //@}

      //@{ \name Initialization methods.

      //! Initialize every element of this vector with \c gamma.
      NOX::Abstract::Vector& init(double gamma);

      //! Initialize every element of this vector with random values
      NOX::Abstract::Vector& random(bool useSeed=false, int seed=1);

      //! Copy source vector \c y into this vector.
      NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector& y);
      NOX::Abstract::Vector& operator=(const NOX::Thyra::Vector& y);

      /*!
       * \brief Put element-wise absolute values of source vector \c y into
       * this vector.
       */
      NOX::Abstract::Vector& abs(const NOX::Abstract::Vector& y);

      //! Put element-wise reciprocal of source vector \c y into this vector.
      NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector& y);

      //@}

      //@{ \name Update methods.

      //! Scale each element of this vector by \c gamma.
      NOX::Abstract::Vector& scale(double gamma);

      //! Scale this vector <em>element-by-element</em> by the vector a.
      NOX::Abstract::Vector& scale(const NOX::Abstract::Vector& a);

      //! Compute x = (alpha * a) + (gamma * x) where x is this vector.
      NOX::Abstract::Vector& update(double alpha,
                    const NOX::Abstract::Vector& a,
                    double gamma = 0.0);

      /*!
       * \brief Compute x = (alpha * a) + (beta * b) + (gamma * x) where x is
       * this vector.
       */
      NOX::Abstract::Vector& update(double alpha,
                    const NOX::Abstract::Vector& a,
                    double beta,
                    const NOX::Abstract::Vector& b,
                    double gamma = 0.0);

      //@}

      //@{ \name Creating new Vectors.

      /*!
       * \brief Create a new %Vector of the same underlying type by
       * cloning "this", and return a pointer to the new vector.
       */
      Teuchos::RCP<NOX::Abstract::Vector>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      /*!
       * \brief Create a MultiVector with \c numVecs+1 columns out of an array
       * of Vectors.  The vector stored under \c this will be the first column
       * with the remaining \c numVecs columns given by \c vecs.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      createMultiVector(const NOX::Abstract::Vector* const* vecs,
            int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

      /*!
       * \brief Create a MultiVector with \c numVecs columns.
       */
      virtual Teuchos::RCP<NOX::Abstract::MultiVector>
      createMultiVector(int numVecs, NOX::CopyType type = NOX::DeepCopy) const;

      //@}

      //@{ \name Norms.

      //!! Norm
      double norm(NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm) const;

      //! Weighted 2-Norm.
      double norm(const NOX::Abstract::Vector& weights) const;

      //@}

      //@{ \name Inner products

      //! Inner product with \c y.
      double innerProduct(const NOX::Abstract::Vector& y) const;

      //@}

      //! Return the length of vector.
      NOX::size_type length() const;

      //! Print the vector.  To be used for debugging only.
      void print(std::ostream& stream) const;

      //! Set the weighting vector used for inner products and norms.
      void setWeightVector(const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weightVec);

      bool hasWeightVector() const;

      //! Returns a weighting vector if one was set, otherwise throws.
      Teuchos::RCP<const ::Thyra::VectorBase<double> >
      getWeightVector() const;

      //@{ \name Derived from NOX::Abstract::ImplicitWeighting

      bool getImplicitWeighting() const;

      void setImplicitWeighting(bool do_implicit_weighting);

      //@}

    protected:

      //! Pointer to Thyra vector owned by this object
      Teuchos::RCP< ::Thyra::VectorBase<double> > thyraVec;

      //! Thyra vector used for weighting inner product and norms
      Teuchos::RCP<const ::Thyra::VectorBase<double> > weightVec_;

      //! Thyra vector used for weighting inner product and norms
      Teuchos::RCP< ::Thyra::VectorBase<double> > tmpVec_;

      //! True if implicit weighting is enabled (i.e. a nonnull wieghtVec_)
      bool do_implicit_weighting_;

    };


  } // namespace Thyra
} // namespace NOX

#endif
