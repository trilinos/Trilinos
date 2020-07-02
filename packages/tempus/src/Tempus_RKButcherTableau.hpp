// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_RKButcherTableau_hpp
#define Tempus_RKButcherTableau_hpp

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include "Tempus_String_Utilities.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


namespace Tempus {


/** \brief Runge-Kutta methods.
 *
 *  This base class specifies the Butcher tableau which defines the
 *  Runge-Kutta (RK) method.  Both explicit and implicit RK methods
 *  can be specified here, and of arbitrary number of stages and orders.
 *  Embedded methods are also supported.
 *
 *  Since this is a generic RK class, no low-storage methods are
 *  incorporated here, however any RK method with a Butcher tableau
 *  can be created with the base class.
 *
 *  There are over 40 derived RK methods that have been implemented,
 *  ranging from first order and eight order, and from single stage
 *  to 5 stages.
 *
 *  This class was taken and modified from Rythmos' RKButcherTableau class.
 */
template<class Scalar>
class RKButcherTableau :
  virtual public Teuchos::Describable,
  virtual public Teuchos::VerboseObject<RKButcherTableau<Scalar> >
{
  public:

    RKButcherTableau(
      std::string stepperType,
      const Teuchos::SerialDenseMatrix<int,Scalar>& A,
      const Teuchos::SerialDenseVector<int,Scalar>& b,
      const Teuchos::SerialDenseVector<int,Scalar>& c,
      const int order,
      const int orderMin,
      const int orderMax,
      const Teuchos::SerialDenseVector<int,Scalar>&
        bstar = Teuchos::SerialDenseVector<int,Scalar>(),
      bool checkC = true)
      : description_(stepperType)
    {
      const int numStages = A.numRows();
      TEUCHOS_ASSERT_EQUALITY( A.numCols(), numStages );
      TEUCHOS_ASSERT_EQUALITY( b.length(), numStages );
      TEUCHOS_ASSERT_EQUALITY( c.length(), numStages );
      TEUCHOS_ASSERT( order > 0 );
      A_ = A;
      b_ = b;
      c_ = c;
      order_ = order;
      orderMin_ = orderMin;
      orderMax_ = orderMax;
      this->set_isImplicit();
      this->set_isDIRK();

      // Consistency check on b
      typedef Teuchos::ScalarTraits<Scalar> ST;
      Scalar sumb = ST::zero();
      for (size_t i = 0; i < this->numStages(); i++) sumb += b_(i);
      TEUCHOS_TEST_FOR_EXCEPTION( std::abs(ST::one()-sumb) > 1.0e-08,
          std::runtime_error,
          "Error - Butcher Tableau b fails to satisfy Sum(b_i) = 1.\n"
          << "          Sum(b_i) = " << sumb << "\n");

      // Consistency check on c
      if (checkC) {
        for (size_t i = 0; i < this->numStages(); i++) {
          Scalar sumai = ST::zero();
          for (size_t j = 0; j < this->numStages(); j++) sumai += A_(i,j);
          bool failed = false;
          if (std::abs(sumai) > 1.0e-08)
            failed = (std::abs((sumai-c_(i))/sumai) > 1.0e-08);
          else
            failed = (std::abs(c_(i)) > 1.0e-08);

          TEUCHOS_TEST_FOR_EXCEPTION( failed, std::runtime_error,
               "Error - Butcher Tableau fails to satisfy c_i = Sum_j(a_ij).\n"
            << "  Stage i      = " << i << "\n"
            << "    c_i         = " << c_(i) << "\n"
            << "    Sum_j(a_ij) = " << sumai << "\n"
            << "  This may be OK as some valid tableaus do not satisfy\n"
            << "  this condition.  If so, construct this RKButcherTableau\n"
            << "  with checkC = false.\n");
        }
      }

      if ( bstar.length() > 0 ) {
        TEUCHOS_ASSERT_EQUALITY( bstar.length(), numStages );
        isEmbedded_ = true;
      } else {
        isEmbedded_ = false;
      }
      bstar_ = bstar;
    }

    /** \brief Return the number of stages */
    virtual std::size_t numStages() const { return A_.numRows(); }
    /** \brief Return the matrix coefficients */
    virtual const Teuchos::SerialDenseMatrix<int,Scalar>& A() const
      { return A_; }
    /** \brief Return the vector of quadrature weights */
    virtual const Teuchos::SerialDenseVector<int,Scalar>& b() const
      { return b_; }
    /** \brief Return the vector of quadrature weights for embedded methods */
    virtual const Teuchos::SerialDenseVector<int,Scalar>& bstar() const
      { return bstar_ ; }
    /** \brief Return the vector of stage positions */
    virtual const Teuchos::SerialDenseVector<int,Scalar>& c() const
      { return c_; }
    /** \brief Return the order */
    virtual int order() const { return order_; }
    /** \brief Return the minimum order */
    virtual int orderMin() const { return orderMin_; }
    /** \brief Return the maximum order */
    virtual int orderMax() const { return orderMax_; }
    /** \brief Return true if the RK method is implicit */
    virtual bool isImplicit() const { return isImplicit_; }
    /** \brief Return true if the RK method is Diagonally Implicit */
    virtual bool isDIRK() const { return isDIRK_; }
    /** \brief Return true if the RK method has embedded capabilities */
    virtual bool isEmbedded() const { return isEmbedded_; }

    /* \brief Redefined from Teuchos::Describable */
    //@{
      virtual std::string description() const { return description_; }

      virtual void describe( Teuchos::FancyOStream &out,
                             const Teuchos::EVerbosityLevel verbLevel) const
      {
        if (verbLevel != Teuchos::VERB_NONE) {
          out << this->description() << std::endl;
          out << "number of Stages = " << this->numStages() << std::endl;
          out << "A = " << printMat(this->A()) << std::endl;
          out << "b = " << printMat(this->b()) << std::endl;
          out << "c = " << printMat(this->c()) << std::endl;
          out << "bstar = " << printMat(this->bstar()) << std::endl;
          out << "order    = " << this->order()    << std::endl;
          out << "orderMin = " << this->orderMin() << std::endl;
          out << "orderMax = " << this->orderMax() << std::endl;
          out << "isImplicit = " << this->isImplicit() << std::endl;
          out << "isDIRK     = " << this->isDIRK()     << std::endl;
          out << "isEmbedded = " << this->isEmbedded() << std::endl;
        }
      }
    //@}

    bool operator == (const RKButcherTableau & t) const {
      const Scalar relTol = 1.0e-15;
      if ( A_->numRows() != t.A_->numRows() ||
           A_->numCols() != t.A_->numCols()    ) {
        return false;
      } else {
        int i, j;
        for(i = 0; i < A_->numRows(); i++) {
          for(j = 0; j < A_->numCols(); j++) {
            if(std::abs((t.A_(i,j) - A_(i,j))/A_(i,j)) > relTol) return false;
          }
        }
      }

      if ( b_->length() != t.b_->length() ||
           b_->length() != t.c_->length() ||
           b_->length() != t.bstar_->length() ) {
        return false;
      } else {
        int i;
        for(i = 0; i < A_->numRows(); i++) {
          if(std::abs((t.b_(i) - b_(i))/b_(i)) > relTol) return false;
          if(std::abs((t.c_(i) - c_(i))/c_(i)) > relTol) return false;
          if(std::abs((t.bstar_(i) - bstar_(i))/bstar_(i)) > relTol) return false;
        }
      }
      return true;
    }

    bool operator != (const RKButcherTableau & t) const {
      return !((*this) == t);
    }

  private:

    RKButcherTableau();

  protected:

    void set_isImplicit() {
      isImplicit_ = false;
      for (size_t i = 0; i < this->numStages(); i++)
        for (size_t j = i; j < this->numStages(); j++)
          if (A_(i,j) != 0.0) isImplicit_ = true;
    }

    /// DIRK is defined as if a_ij = 0 for j>i and a_ii != 0 for at least one i.
    void set_isDIRK() {
      isDIRK_ = true;
      bool nonZero = false;
      for (size_t i = 0; i < this->numStages(); i++) {
        if (A_(i,i) != 0.0) nonZero = true;
        for (size_t j = i+1; j < this->numStages(); j++)
          if (A_(i,j) != 0.0) isDIRK_ = false;
      }
      if (nonZero == false) isDIRK_ = false;
    }

    std::string description_;

    Teuchos::SerialDenseMatrix<int,Scalar> A_;
    Teuchos::SerialDenseVector<int,Scalar> b_;
    Teuchos::SerialDenseVector<int,Scalar> c_;
    int order_;
    int orderMin_;
    int orderMax_;
    bool isImplicit_;
    bool isDIRK_;
    bool isEmbedded_;
    Teuchos::SerialDenseVector<int,Scalar> bstar_;
};


} // namespace Tempus


#endif // Tempus_RKButcherTableau_hpp
