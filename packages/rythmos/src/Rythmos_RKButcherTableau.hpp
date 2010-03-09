//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER


#ifndef RYTHMOS_RK_BUTCHER_TABLEAU_HPP
#define RYTHMOS_RK_BUTCHER_TABLEAU_HPP

#include "Rythmos_Types.hpp"
#include "Rythmos_RKButcherTableauBase.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

#include "Thyra_ProductVectorBase.hpp"

namespace Rythmos {

  using Teuchos::as;

  inline const std::string RKBT_ForwardEuler_name() { return  "Forward Euler"; } // done
  inline const std::string RKBT_BackwardEuler_name() { return  "Backward Euler"; } // done
  inline const std::string Explicit4Stage_name() { return  "Explicit 4 Stage"; } // done
  inline const std::string Explicit3_8Rule_name() { return  "Explicit 3/8 Rule"; } // done

  inline const std::string Explicit2Stage2ndOrderRunge_name() { return  "Explicit 2 Stage 2nd order by Runge"; } // done
  inline const std::string Explicit3Stage3rdOrderHeun_name() { return  "Explicit 3 Stage 3rd order by Heun"; } // done
  inline const std::string Explicit3Stage3rdOrder_name() { return  "Explicit 3 Stage 3rd order"; } // done
  inline const std::string Explicit4Stage3rdOrderRunge_name() { return  "Explicit 4 Stage 3rd order by Runge"; } // done

  inline const std::string Implicit1Stage2ndOrderGauss_name() { return  "Implicit 1 Stage 2nd order Gauss"; } // done
  inline const std::string Implicit2Stage4thOrderGauss_name() { return  "Implicit 2 Stage 4th order Gauss"; } // done
  inline const std::string Implicit3Stage6thOrderGauss_name() { return  "Implicit 3 Stage 6th order Gauss"; } // done

  inline const std::string Implicit1Stage1stOrderRadauA_name() { return  "Implicit 1 Stage 1st order Radau left"; } // done
  inline const std::string Implicit2Stage3rdOrderRadauA_name() { return  "Implicit 2 Stage 3rd order Radau left"; } // done
  inline const std::string Implicit3Stage5thOrderRadauA_name() { return  "Implicit 3 Stage 5th order Radau left"; } // done

  inline const std::string Implicit1Stage1stOrderRadauB_name() { return  "Implicit 1 Stage 1st order Radau right"; } // done
  inline const std::string Implicit2Stage3rdOrderRadauB_name() { return  "Implicit 2 Stage 3rd order Radau right"; } // done
  inline const std::string Implicit3Stage5thOrderRadauB_name() { return  "Implicit 3 Stage 5th order Radau right"; } // done

  inline const std::string Implicit2Stage2ndOrderLobattoA_name() { return  "Implicit 2 Stage 2nd order Lobatto A"; } // done
  inline const std::string Implicit3Stage4thOrderLobattoA_name() { return  "Implicit 3 Stage 4th order Lobatto A"; } // done
  inline const std::string Implicit4Stage6thOrderLobattoA_name() { return  "Implicit 4 Stage 6th order Lobatto A"; } // done

  inline const std::string Implicit2Stage2ndOrderLobattoB_name() { return  "Implicit 2 Stage 2nd order Lobatto B"; } // done
  inline const std::string Implicit3Stage4thOrderLobattoB_name() { return  "Implicit 3 Stage 4th order Lobatto B"; } // done
  inline const std::string Implicit4Stage6thOrderLobattoB_name() { return  "Implicit 4 Stage 6th order Lobatto B"; } // done

  inline const std::string Implicit2Stage2ndOrderLobattoC_name() { return  "Implicit 2 Stage 2nd order Lobatto C"; } // done
  inline const std::string Implicit3Stage4thOrderLobattoC_name() { return  "Implicit 3 Stage 4th order Lobatto C"; } // done
  inline const std::string Implicit4Stage6thOrderLobattoC_name() { return  "Implicit 4 Stage 6th order Lobatto C"; } // done

  inline const std::string Implicit2Stage4thOrderHammerHollingsworth_name() { return  "Implicit 2 Stage 4th Order Hammer & Hollingsworth"; } // done
  inline const std::string Implicit3Stage6thOrderKuntzmannButcher_name() { return  "Implicit 3 Stage 6th Order Kuntzmann & Butcher"; } // done
  inline const std::string Implicit4Stage8thOrderKuntzmannButcher_name() { return  "Implicit 4 Stage 8th Order Kuntzmann & Butcher"; } // done

  inline const std::string DIRK2Stage3rdOrder_name() { return  "Diagonal IRK 2 Stage 3rd order"; } // done

  inline const std::string SDIRK2Stage3rdOrder_name() { return  "Singly Diagonal IRK 2 Stage 3rd order"; } // done
  inline const std::string SDIRK5Stage5thOrder_name() { return  "Singly Diagonal IRK 5 Stage 5th order"; } // done
  inline const std::string SDIRK5Stage4thOrder_name() { return  "Singly Diagonal IRK 5 Stage 4th order"; } // done
  inline const std::string SDIRK3Stage4thOrder_name() { return  "Singly Diagonal IRK 3 Stage 4th order"; } // done

template<class Scalar>
class RKButcherTableauDefaultBase :
  virtual public RKButcherTableauBase<Scalar>,
  virtual public Teuchos::ParameterListAcceptorDefaultBase
{
  public:
    /** \brief . */
    virtual int numStages() const { return A_.numRows(); }
    /** \brief . */
    virtual const Teuchos::SerialDenseMatrix<int,Scalar>& A() const { return A_; }
    /** \brief . */
    virtual const Teuchos::SerialDenseVector<int,Scalar>& b() const { return b_; }
    /** \brief . */
    virtual const Teuchos::SerialDenseVector<int,Scalar>& c() const { return c_; }
    /** \brief . */
    virtual int order() const { return order_; }
    /** \brief . */
    virtual void setDescription(std::string longDescription) { longDescription_ = longDescription; }

    /** \brief . */
    virtual void initialize(
        const Teuchos::SerialDenseMatrix<int,Scalar>& A_in,
        const Teuchos::SerialDenseVector<int,Scalar>& b_in,
        const Teuchos::SerialDenseVector<int,Scalar>& c_in,
        const int order_in,
        const std::string& longDescription_in
        )
    { 
      const int numStages_in = A_in.numRows();
      TEUCHOS_ASSERT_EQUALITY( A_in.numRows(), numStages_in );
      TEUCHOS_ASSERT_EQUALITY( A_in.numCols(), numStages_in );
      TEUCHOS_ASSERT_EQUALITY( b_in.length(), numStages_in );
      TEUCHOS_ASSERT_EQUALITY( c_in.length(), numStages_in );
      TEUCHOS_ASSERT( order_in > 0 );
      A_ = A_in;
      b_ = b_in;
      c_ = c_in;
      order_ = order_in;
      longDescription_ = longDescription_in;
    }

    /* \brief Redefined from Teuchos::ParameterListAcceptorDefaultBase */
    //@{
    
    /** \brief . */
    virtual void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
    {
      TEST_FOR_EXCEPT( is_null(paramList) );
      paramList->validateParameters(*this->getValidParameters());
      Teuchos::readVerboseObjectSublist(&*paramList,this);
      setMyParamList(paramList);
    }

    /** \brief . */
    virtual RCP<const Teuchos::ParameterList> getValidParameters() const 
    { 
      if (is_null(validPL_)) {
        validPL_ = Teuchos::parameterList();
        validPL_->set("Description","",this->getMyDescription());
        Teuchos::setupVerboseObjectSublist(&*validPL_);
      }
      return validPL_; 
    }

    //@}

    /* \brief Redefined from Teuchos::Describable */
    //@{
    
    /** \brief . */
    virtual std::string description() const { return "Rythmos::RKButcherTableauDefaultBase"; }

    /** \brief . */
    virtual void describe(
      Teuchos::FancyOStream &out,
      const Teuchos::EVerbosityLevel verbLevel
      ) const
    {
      if (verbLevel != Teuchos::VERB_NONE) {
        out << this->description() << std::endl;
        out << this->getMyDescription() << std::endl;
        out << "number of Stages = " << this->numStages() << std::endl;
        out << "A = " << this->A() << std::endl;
        out << "b = " << this->b() << std::endl;
        out << "c = " << this->c() << std::endl;
        out << "order = " << this->order() << std::endl;
      }
    }

    //@}

  protected:
    void setMyDescription(std::string longDescription) { longDescription_ = longDescription; }
    const std::string& getMyDescription() const { return longDescription_; }

    void setMy_A(const Teuchos::SerialDenseMatrix<int,Scalar>& new_A) { A_ = new_A; }
    void setMy_b(const Teuchos::SerialDenseVector<int,Scalar>& new_b) { b_ = new_b; }
    void setMy_c(const Teuchos::SerialDenseVector<int,Scalar>& new_c) { c_ = new_c; }
    void setMy_order(const int& new_order) { order_ = new_order; }

    void setMyValidParameterList( const RCP<ParameterList> validPL ) { validPL_ = validPL; }
    RCP<ParameterList> getMyNonconstValidParameterList() { return validPL_; }

  private:
    Teuchos::SerialDenseMatrix<int,Scalar> A_;
    Teuchos::SerialDenseVector<int,Scalar> b_;
    Teuchos::SerialDenseVector<int,Scalar> c_;
    int order_;
    std::string longDescription_;
    mutable RCP<ParameterList> validPL_;
};


// Nonmember constructor
template<class Scalar>
RCP<RKButcherTableauBase<Scalar> > rKButcherTableau()
{
  return(rcp(new RKButcherTableauDefaultBase<Scalar>()));
}

// Nonmember constructor
template<class Scalar>
RCP<RKButcherTableauBase<Scalar> > rKButcherTableau(
    const Teuchos::SerialDenseMatrix<int,Scalar>& A_in,
    const Teuchos::SerialDenseVector<int,Scalar>& b_in,
    const Teuchos::SerialDenseVector<int,Scalar>& c_in,
    int order_in,
    const std::string& description_in = ""
    )
{
  RCP<RKButcherTableauDefaultBase<Scalar> > rkbt = rcp(new RKButcherTableauDefaultBase<Scalar>());
  rkbt->initialize(A_in,b_in,c_in,order_in,description_in);
  return(rkbt);
}


template<class Scalar>
class BackwardEuler_RKBT : 
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
  BackwardEuler_RKBT() 
  {
    std::ostringstream description;
    description << RKBT_BackwardEuler_name() << "\n"
                << "c = [ 1 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'" << std::endl;
    typedef ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    A(0,0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    b(0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    c(0) = ST::one();

    this->setMyDescription(description.str());
    this->setMy_A(A);
    this->setMy_b(b);
    this->setMy_c(c);
    this->setMy_order(1);
  }
};



template<class Scalar>
class ForwardEuler_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:

    ForwardEuler_RKBT()
    {
      std::ostringstream description;
      description << RKBT_ForwardEuler_name() << "\n"
                  << "c = [ 0 ]'\n"
                  << "A = [ 0 ]\n"
                  << "b = [ 1 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
      Teuchos::SerialDenseVector<int,Scalar> b(1);
      b(0) = ST::one();
      Teuchos::SerialDenseVector<int,Scalar> c(1);

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(1);
    }
};


template<class Scalar>
class Explicit4Stage4thOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit4Stage4thOrder_RKBT() 
    {
      std::ostringstream description;
      description << Explicit4Stage_name() << "\n"
                  << "\"The\" Runge-Kutta Method (explicit):\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Edition\n"
                  << "E. Hairer, S.P. Norsett, G. Wanner\n"
                  << "Table 1.2, pg 138\n"
                  << "c = [  0  1/2 1/2  1  ]'\n"
                  << "A = [  0              ] \n"
                  << "    [ 1/2  0          ]\n"
                  << "    [  0  1/2  0      ]\n"
                  << "    [  0   0   1   0  ]\n"
                  << "b = [ 1/6 1/3 1/3 1/6 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      Scalar onehalf = ST::one()/(2*ST::one());
      Scalar onesixth = ST::one()/(6*ST::one());
      Scalar onethird = ST::one()/(3*ST::one());

      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;

      A(1,0) = onehalf;
      A(1,1) = zero;
      A(1,2) = zero;
      A(1,3) = zero;

      A(2,0) = zero;
      A(2,1) = onehalf;
      A(2,2) = zero;
      A(2,3) = zero;

      A(3,0) = zero;
      A(3,1) = zero;
      A(3,2) = one;
      A(3,3) = zero;

      // Fill b:
      b(0) = onesixth;
      b(1) = onethird;
      b(2) = onethird;
      b(3) = onesixth;
      
      // fill b_c_
      c(0) = zero;
      c(1) = onehalf;
      c(2) = onehalf;
      c(3) = one;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Explicit3_8Rule_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit3_8Rule_RKBT()
    { 

      std::ostringstream description;
      description << Explicit3_8Rule_name() << "\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Edition\n"
                  << "E. Hairer, S.P. Norsett, G. Wanner\n"
                  << "Table 1.2, pg 138\n"
                  << "c = [  0  1/3 2/3  1  ]'\n"
                  << "A = [  0              ]\n"
                  << "    [ 1/3  0          ]\n"
                  << "    [-1/3  1   0      ]\n"
                  << "    [  1  -1   1   0  ]\n"
                  << "b = [ 1/8 3/8 3/8 1/8 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      Scalar one = ST::one();
      Scalar zero = ST::zero();
      Scalar one_third    = as<Scalar>(ST::one()/(3*ST::one()));
      Scalar two_third    = as<Scalar>(2*ST::one()/(3*ST::one()));
      Scalar one_eighth   = as<Scalar>(ST::one()/(8*ST::one()));
      Scalar three_eighth = as<Scalar>(3*ST::one()/(8*ST::one()));

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;

      A(1,0) = one_third;
      A(1,1) = zero;
      A(1,2) = zero;
      A(1,3) = zero;

      A(2,0) = as<Scalar>(-one_third);
      A(2,1) = one;
      A(2,2) = zero;
      A(2,3) = zero;

      A(3,0) = one;
      A(3,1) = as<Scalar>(-one);
      A(3,2) = one;
      A(3,3) = zero;

      // Fill b:
      b(0) = one_eighth;
      b(1) = three_eighth;
      b(2) = three_eighth;
      b(3) = one_eighth;
      
      // Fill c:
      c(0) = zero;
      c(1) = one_third;
      c(2) = two_third;
      c(3) = one;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Explicit4Stage3rdOrderRunge_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit4Stage3rdOrderRunge_RKBT()
    { 

      std::ostringstream description;
      description << Explicit4Stage3rdOrderRunge_name() << "\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Edition\n"
                  << "E. Hairer, S.P. Norsett, G. Wanner\n"
                  << "Table 1.1, pg 135\n"
                  << "c = [  0  1/2  1   1  ]'\n"
                  << "A = [  0              ]\n"
                  << "    [ 1/2  0          ]\n"
                  << "    [  0   1   0      ]\n"
                  << "    [  0   0   1   0  ]\n"
                  << "b = [ 1/6 2/3  0  1/6 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      Scalar one = ST::one();
      Scalar onehalf = ST::one()/(2*ST::one());
      Scalar onesixth = ST::one()/(6*ST::one());
      Scalar twothirds = 2*ST::one()/(3*ST::one());
      Scalar zero = ST::zero();

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;

      A(1,0) = onehalf;
      A(1,1) = zero;
      A(1,2) = zero;
      A(1,3) = zero;

      A(2,0) = zero;
      A(2,1) = one;
      A(2,2) = zero;
      A(2,3) = zero;

      A(3,0) = zero;
      A(3,1) = zero;
      A(3,2) = one;
      A(3,3) = zero;

      // Fill b:
      b(0) = onesixth;
      b(1) = twothirds;
      b(2) = zero;
      b(3) = onesixth;
      
      // Fill c:
      c(0) = zero;
      c(1) = onehalf;
      c(2) = one;
      c(3) = one;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Explicit3Stage3rdOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit3Stage3rdOrder_RKBT()
    { 

      std::ostringstream description;
      description << Explicit3Stage3rdOrder_name() << "\n" 
                  << "c = [  0  1/2  1  ]'\n"
                  << "A = [  0          ]\n"
                  << "    [ 1/2  0      ]\n"
                  << "    [ -1   2   0  ]\n"
                  << "b = [ 1/6 4/6 1/6 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      Scalar one = ST::one();
      Scalar two = as<Scalar>(2*ST::one());
      Scalar zero = ST::zero();
      Scalar onehalf = ST::one()/(2*ST::one());
      Scalar onesixth = ST::one()/(6*ST::one());
      Scalar foursixth = 4*ST::one()/(6*ST::one());

      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;

      A(1,0) = onehalf;
      A(1,1) = zero;
      A(1,2) = zero;

      A(2,0) = -one;
      A(2,1) = two;
      A(2,2) = zero;

      // Fill b:
      b(0) = onesixth;
      b(1) = foursixth;
      b(2) = onesixth;
      
      // fill b_c_
      c(0) = zero;
      c(1) = onehalf;
      c(2) = one;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Explicit3Stage3rdOrderHeun_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit3Stage3rdOrderHeun_RKBT()
    { 
      std::ostringstream description;
      description << Explicit3Stage3rdOrderHeun_name() << "\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Edition\n"
                  << "E. Hairer, S.P. Norsett, G. Wanner\n"
                  << "Table 1.1, pg 135\n"
                  << "c = [  0  1/3 2/3 ]'\n"
                  << "A = [  0          ] \n"
                  << "    [ 1/3  0      ]\n"
                  << "    [  0  2/3  0  ]\n"
                  << "b = [ 1/4  0  3/4 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      Scalar onethird = one/(3*one);
      Scalar twothirds = 2*one/(3*one);
      Scalar onefourth = one/(4*one);
      Scalar threefourths = 3*one/(4*one);

      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;

      A(1,0) = onethird;
      A(1,1) = zero;
      A(1,2) = zero;

      A(2,0) = zero;
      A(2,1) = twothirds;
      A(2,2) = zero;

      // Fill b:
      b(0) = onefourth;
      b(1) = zero;
      b(2) = threefourths;
      
      // fill b_c_
      c(0) = zero;
      c(1) = onethird;
      c(2) = twothirds;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Explicit2Stage2ndOrderRunge_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Explicit2Stage2ndOrderRunge_RKBT()
    { 
      std::ostringstream description;
      description << Explicit2Stage2ndOrderRunge_name() << "\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Edition\n"
                  << "E. Hairer, S.P. Norsett, G. Wanner\n"
                  << "Table 1.1, pg 135\n"
                  << "c = [  0  1/2 ]'\n"
                  << "A = [  0      ]\n"
                  << "    [ 1/2  0  ]\n"
                  << "b = [  0   1  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      Scalar onehalf = ST::one()/(2*ST::one());

      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);

      // Fill A:
      A(0,0) = zero;
      A(0,1) = zero;

      A(1,0) = onehalf;
      A(1,1) = zero;

      // Fill b:
      b(0) = zero;
      b(1) = one;
      
      // fill b_c_
      c(0) = zero;
      c(1) = onehalf;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(2);
    }
};


// 04/07/09 tscoffe:  I verified manually that the Convergence Testing passes
// with gamma_default_ = -1.
template<class Scalar>
class SDIRK2Stage3rdOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    SDIRK2Stage3rdOrder_RKBT()
    { 
      std::ostringstream description;
      description << SDIRK2Stage3rdOrder_name() << "\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.2, pg 207\n"
                  << "gamma = (3+-sqrt(3))/6\n"
                  << "c = [  gamma     1-gamma  ]'\n"
                  << "A = [  gamma     0        ]\n"
                  << "    [ 1-2*gamma  gamma    ]\n"
                  << "b = [ 1/2        1/2      ]'" << std::endl;

      this->setMyDescription(description.str());
      gamma_coeff_default_ = 1;
      gamma_coeff_ = gamma_coeff_default_;
      this->setupData();

      RCP<ParameterList> validPL = Teuchos::parameterList();
      validPL->set("Description","",this->getMyDescription());
      validPL->set<int>("gamma coefficient",gamma_coeff_default_,"gamma = (3+[gamma coefficient]*sqrt(3))/6, [gamma coefficient] = +-1");
      Teuchos::setupVerboseObjectSublist(&*validPL);
      this->setMyValidParameterList(validPL);
    }
    void setupData() 
    {
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      Scalar gamma = as<Scalar>( (3*one + as<Scalar>(gamma_coeff_)*ST::squareroot(3*one))/(6*one) );
      A(0,0) = gamma;
      A(0,1) = zero;
      A(1,0) = as<Scalar>( one - 2*gamma );
      A(1,1) = gamma;
      b(0) = as<Scalar>( one/(2*one) );
      b(1) = as<Scalar>( one/(2*one) );
      c(0) = gamma;
      c(1) = as<Scalar>( one - gamma );

      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
    void setParameterList(RCP<Teuchos::ParameterList> const& paramList)
    {
      TEST_FOR_EXCEPT( is_null(paramList) );
      paramList->validateParameters(*this->getValidParameters());
      Teuchos::readVerboseObjectSublist(&*paramList,this);
      gamma_coeff_ = paramList->get<int>("gamma coefficient",gamma_coeff_default_);
      this->setupData();
      this->setMyParamList(paramList);
    }
  private:
    int gamma_coeff_default_;
    int gamma_coeff_;
};


template<class Scalar>
class DIRK2Stage3rdOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    DIRK2Stage3rdOrder_RKBT()
    { 

      std::ostringstream description;
      description << DIRK2Stage3rdOrder_name() << "\n"
                  << "Hammer & Hollingsworth method\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.1, pg 205\n"
                  << "c = [  0    0  ]'\n"
                  << "A = [  0    0  ]\n"
                  << "    [ 1/3  1/3 ]\n"
                  << "b = [ 1/4  3/4 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      A(0,0) = zero;
      A(0,1) = zero;
      A(1,0) = as<Scalar>( one/(3*one) );
      A(1,1) = as<Scalar>( one/(3*one) );
      b(0) = as<Scalar>( one/(4*one) );
      b(1) = as<Scalar>( 3*one/(4*one) );
      c(0) = zero;
      c(1) = as<Scalar>( 2*one/(3*one) );
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Implicit3Stage6thOrderKuntzmannButcher_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage6thOrderKuntzmannButcher_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage6thOrderKuntzmannButcher_name() << "\n"
                  << "Kuntzmann & Butcher method\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.4, pg 209\n"
                  << "c = [ 1/2-sqrt(15)/10   1/2              1/2-sqrt(15)/10  ]'\n"
                  << "A = [ 5/36              2/9-sqrt(15)/15  5/36-sqrt(15)/30 ]\n"
                  << "    [ 5/36+sqrt(15)/24  2/9              5/36-sqrt(15)/24 ]\n"
                  << "    [ 5/36+sqrt(15)/30  2/9+sqrt(15)/15  5/36             ]\n"
                  << "b = [ 5/18              4/9              5/18             ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( 5*one/(36*one) );
      A(0,1) = as<Scalar>( 2*one/(9*one) - ST::squareroot(15*one)/(15*one) );
      A(0,2) = as<Scalar>( 5*one/(36*one) - ST::squareroot(15*one)/(30*one) );
      A(1,0) = as<Scalar>( 5*one/(36*one) + ST::squareroot(15*one)/(24*one) );
      A(1,1) = as<Scalar>( 2*one/(9*one) );
      A(1,2) = as<Scalar>( 5*one/(36*one) - ST::squareroot(15*one)/(24*one) );
      A(2,0) = as<Scalar>( 5*one/(36*one) + ST::squareroot(15*one)/(30*one) );
      A(2,1) = as<Scalar>( 2*one/(9*one) + ST::squareroot(15*one)/(15*one) );
      A(2,2) = as<Scalar>( 5*one/(36*one) );
      b(0) = as<Scalar>( 5*one/(18*one) );
      b(1) = as<Scalar>( 4*one/(9*one) );
      b(2) = as<Scalar>( 5*one/(18*one) );
      c(0) = as<Scalar>( one/(2*one)-ST::squareroot(15*one)/(10*one) );
      c(1) = as<Scalar>( one/(2*one) );
      c(2) = as<Scalar>( one/(2*one)+ST::squareroot(15*one)/(10*one) );
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(6);
    }
};


template<class Scalar>
class Implicit4Stage8thOrderKuntzmannButcher_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit4Stage8thOrderKuntzmannButcher_RKBT()
    { 

      std::ostringstream description;
      description << Implicit4Stage8thOrderKuntzmannButcher_name() << "\n"
                  << "Kuntzmann & Butcher method\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.5, pg 209\n"
                  << "c = [ 1/2-w2     1/2-w2p     1/2+w2p     1/2+w2    ]'\n"
                  << "A = [ w1         w1p-w3+w4p  w1p-w3-w4p  w1-w5     ]\n"
                  << "    [ w1-w3p+w4  w1p         w1p-w5p     w1-w3p-w4 ]\n"
                  << "    [ w1+w3p+w4  w1p+w5p     w1p         w1+w3p-w4 ]\n"
                  << "    [ w1+w5      w1p+w3+w4p  w1p+w3-w4p  w1        ]\n"
                  << "b = [ 2*w1       2*w1p       2*w1p       2*w1      ]'\n"
                  << "w1 = 1/8-sqrt(30)/144\n"
                  << "w2 = (1/2)*sqrt((15+2*sqrt(3))/35)\n"
                  << "w3 = w2*(1/6+sqrt(30)/24)\n"
                  << "w4 = w2*(1/21+5*sqrt(30)/168)\n"
                  << "w5 = w2-2*w3\n"
                  << "w1p = 1/8+sqrt(30)/144\n"
                  << "w2p = (1/2)*sqrt((15-2*sqrt(3))/35)\n"
                  << "w3p = w2*(1/6-sqrt(30)/24)\n"
                  << "w4p = w2*(1/21-5*sqrt(30)/168)\n"
                  << "w5p = w2p-2*w3p" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar onehalf = as<Scalar>( one/(2*one) );
      Scalar w1 = as<Scalar>( one/(8*one) - ST::squareroot(30*one)/(144*one) );
      Scalar w2 = as<Scalar>( (one/(2*one))*ST::squareroot((15*one+2*one*ST::squareroot(30*one))/(35*one)) );
      Scalar w3 = as<Scalar>( w2*(one/(6*one)+ST::squareroot(30*one)/(24*one)) );
      Scalar w4 = as<Scalar>( w2*(one/(21*one)+5*one*ST::squareroot(30*one)/(168*one)) );
      Scalar w5 = as<Scalar>( w2-2*w3 );
      Scalar w1p = as<Scalar>( one/(8*one) + ST::squareroot(30*one)/(144*one) );
      Scalar w2p = as<Scalar>( (one/(2*one))*ST::squareroot((15*one-2*one*ST::squareroot(30*one))/(35*one)) );
      Scalar w3p = as<Scalar>( w2p*(one/(6*one)-ST::squareroot(30*one)/(24*one)) );
      Scalar w4p = as<Scalar>( w2p*(one/(21*one)-5*one*ST::squareroot(30*one)/(168*one)) );
      Scalar w5p = as<Scalar>( w2p-2*w3p );
      A(0,0) = w1;
      A(0,1) = w1p-w3+w4p;
      A(0,2) = w1p-w3-w4p;
      A(0,3) = w1-w5;
      A(1,0) = w1-w3p+w4;
      A(1,1) = w1p;
      A(1,2) = w1p-w5p;
      A(1,3) = w1-w3p-w4;
      A(2,0) = w1+w3p+w4;
      A(2,1) = w1p+w5p;
      A(2,2) = w1p;
      A(2,3) = w1+w3p-w4;
      A(3,0) = w1+w5;
      A(3,1) = w1p+w3+w4p;
      A(3,2) = w1p+w3-w4p;
      A(3,3) = w1;
      b(0) = 2*w1;
      b(1) = 2*w1p;
      b(2) = 2*w1p;
      b(3) = 2*w1; 
      c(0) = onehalf - w2;
      c(1) = onehalf - w2p;
      c(2) = onehalf + w2p;
      c(3) = onehalf + w2;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(8);
    }
};


template<class Scalar>
class Implicit2Stage4thOrderHammerHollingsworth_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage4thOrderHammerHollingsworth_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage4thOrderHammerHollingsworth_name() << "\n"
                  << "Hammer & Hollingsworth method\n"
                  << "Solving Ordinary Differential Equations I:  Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.3, pg 207\n"
                  << "c = [ 1/2-sqrt(3)/6  1/2+sqrt(3)/6 ]'\n"
                  << "A = [ 1/4            1/4-sqrt(3)/6 ]\n"
                  << "    [ 1/4+sqrt(3)/6  1/4           ]\n"
                  << "b = [ 1/2            1/2           ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar onequarter = as<Scalar>( one/(4*one) );
      Scalar onehalf = as<Scalar>( one/(2*one) );
      A(0,0) = onequarter;
      A(0,1) = as<Scalar>( onequarter-ST::squareroot(3*one)/(6*one) );
      A(1,0) = as<Scalar>( onequarter+ST::squareroot(3*one)/(6*one) );
      A(1,1) = onequarter;
      b(0) = onehalf;
      b(1) = onehalf;
      c(0) = as<Scalar>( onehalf - ST::squareroot(3*one)/(6*one) );
      c(1) = as<Scalar>( onehalf + ST::squareroot(3*one)/(6*one) );
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Implicit1Stage2ndOrderGauss_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit1Stage2ndOrderGauss_RKBT()
    { 

      std::ostringstream description;
      description << Implicit1Stage2ndOrderGauss_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.2, pg 72\n"
                  << "Also:  Implicit midpoint rule\n"
                  << "Solving Ordinary Differential Equations I: Nonstiff Problems, 2nd Revised Edition\n"
                  << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                  << "Table 7.1, pg 205\n"
                  << "c = [ 1/2 ]'\n"
                  << "A = [ 1/2 ]\n"
                  << "b = [  1  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 1;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar onehalf = ST::one()/(2*ST::one());
      Scalar one = ST::one();
      A(0,0) = onehalf;
      b(0) = one;
      c(0) = onehalf;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(2);
    }
};


template<class Scalar>
class Implicit2Stage4thOrderGauss_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage4thOrderGauss_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage4thOrderGauss_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.2, pg 72\n"
                  << "c = [ 1/2-sqrt(3)/6  1/2+sqrt(3)/6 ]'\n"
                  << "A = [ 1/4            1/4-sqrt(3)/6 ]\n"
                  << "    [ 1/4+sqrt(3)/6  1/4           ]\n"
                  << "b = [ 1/2            1/2 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar onehalf = as<Scalar>(one/(2*one));
      Scalar three = as<Scalar>(3*one);
      Scalar six = as<Scalar>(6*one);
      Scalar onefourth = as<Scalar>(one/(4*one));
      Scalar alpha = ST::squareroot(three)/six;

      A(0,0) = onefourth;
      A(0,1) = onefourth-alpha;
      A(1,0) = onefourth+alpha;
      A(1,1) = onefourth;
      b(0) = onehalf;
      b(1) = onehalf;
      c(0) = onehalf-alpha;
      c(1) = onehalf+alpha;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Implicit3Stage6thOrderGauss_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage6thOrderGauss_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage6thOrderGauss_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.2, pg 72\n"
                  << "c = [ 1/2-sqrt(15)/10   1/2              1/2+sqrt(15)/10  ]'\n"
                  << "A = [ 5/36              2/9-sqrt(15)/15  5/36-sqrt(15)/30 ]\n"
                  << "    [ 5/36+sqrt(15)/24  2/9              5/36-sqrt(15)/24 ]\n"
                  << "    [ 5/36+sqrt(15)/30  2/9+sqrt(15)/15  5/36             ]\n"
                  << "b = [ 5/18              4/9              5/18             ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar ten = as<Scalar>(10*one);
      Scalar fifteen = as<Scalar>(15*one);
      Scalar twentyfour = as<Scalar>(24*one);
      Scalar thirty = as<Scalar>(30*one);
      Scalar sqrt15over10 = as<Scalar>(ST::squareroot(fifteen)/ten);
      Scalar sqrt15over15 = as<Scalar>(ST::squareroot(fifteen)/fifteen);
      Scalar sqrt15over24 = as<Scalar>(ST::squareroot(fifteen)/twentyfour);
      Scalar sqrt15over30 = as<Scalar>(ST::squareroot(fifteen)/thirty);

      A(0,0) = as<Scalar>(5*one/(36*one));
      A(0,1) = as<Scalar>(2*one/(9*one))-sqrt15over15;
      A(0,2) = as<Scalar>(5*one/(36*one))-sqrt15over30;
      A(1,0) = as<Scalar>(5*one/(36*one))+sqrt15over24;
      A(1,1) = as<Scalar>(2*one/(9*one));
      A(1,2) = as<Scalar>(5*one/(36*one))-sqrt15over24;
      A(2,0) = as<Scalar>(5*one/(36*one))+sqrt15over30;
      A(2,1) = as<Scalar>(2*one/(9*one))+sqrt15over15;
      A(2,2) = as<Scalar>(5*one/(36*one));
      b(0) = as<Scalar>(5*one/(18*one));
      b(1) = as<Scalar>(4*one/(9*one));
      b(2) = as<Scalar>(5*one/(18*one));
      c(0) = as<Scalar>(one/(2*one))-sqrt15over10;
      c(1) = as<Scalar>(one/(2*one));
      c(2) = as<Scalar>(one/(2*one))+sqrt15over10;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(6);
    }
};


template<class Scalar>
class Implicit1Stage1stOrderRadauA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit1Stage1stOrderRadauA_RKBT()
    { 

      std::ostringstream description;
      description << Implicit1Stage1stOrderRadauA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.3, pg 73\n"
                  << "c = [ 0 ]'\n"
                  << "A = [ 1 ]\n"
                  << "b = [ 1 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 1;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      Scalar zero = ST::zero();
      A(0,0) = one;
      b(0) = one;
      c(0) = zero;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(1);
    }
};


template<class Scalar>
class Implicit2Stage3rdOrderRadauA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage3rdOrderRadauA_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage3rdOrderRadauA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.3, pg 73\n"
                  << "c = [  0    2/3 ]'\n"
                  << "A = [ 1/4  -1/4 ]\n"
                  << "    [ 1/4  5/12 ]\n"
                  << "b = [ 1/4  3/4  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>(one/(4*one));
      A(0,1) = as<Scalar>(-one/(4*one));
      A(1,0) = as<Scalar>(one/(4*one));
      A(1,1) = as<Scalar>(5*one/(12*one));
      b(0) = as<Scalar>(one/(4*one));
      b(1) = as<Scalar>(3*one/(4*one));
      c(0) = zero;
      c(1) = as<Scalar>(2*one/(3*one));
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Implicit3Stage5thOrderRadauA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage5thOrderRadauA_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage5thOrderRadauA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.4, pg 73\n"
                  << "c = [  0   (6-sqrt(6))/10       (6+sqrt(6))/10      ]'\n"
                  << "A = [ 1/9  (-1-sqrt(6))/18      (-1+sqrt(6))/18     ]\n"
                  << "    [ 1/9  (88+7*sqrt(6))/360   (88-43*sqrt(6))/360 ]\n"
                  << "    [ 1/9  (88+43*sqrt(6))/360  (88-7*sqrt(6))/360  ]\n"
                  << "b = [ 1/9  (16+sqrt(6))/36      (16-sqrt(6))/36     ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>(one/(9*one));
      A(0,1) = as<Scalar>( (-one-ST::squareroot(6*one))/(18*one) );
      A(0,2) = as<Scalar>( (-one+ST::squareroot(6*one))/(18*one) );
      A(1,0) = as<Scalar>(one/(9*one));
      A(1,1) = as<Scalar>( (88*one+7*one*ST::squareroot(6*one))/(360*one) );
      A(1,2) = as<Scalar>( (88*one-43*one*ST::squareroot(6*one))/(360*one) );
      A(2,0) = as<Scalar>(one/(9*one));
      A(2,1) = as<Scalar>( (88*one+43*one*ST::squareroot(6*one))/(360*one) );
      A(2,2) = as<Scalar>( (88*one-7*one*ST::squareroot(6*one))/(360*one) );
      b(0) = as<Scalar>(one/(9*one));
      b(1) = as<Scalar>( (16*one+ST::squareroot(6*one))/(36*one) );
      b(2) = as<Scalar>( (16*one-ST::squareroot(6*one))/(36*one) );
      c(0) = zero;
      c(1) = as<Scalar>( (6*one-ST::squareroot(6*one))/(10*one) );
      c(2) = as<Scalar>( (6*one+ST::squareroot(6*one))/(10*one) );
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(5);
    }
};


template<class Scalar>
class Implicit1Stage1stOrderRadauB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit1Stage1stOrderRadauB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit1Stage1stOrderRadauB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.5, pg 74\n"
                  << "c = [ 1 ]'\n"
                  << "A = [ 1 ]\n"
                  << "b = [ 1 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 1;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      A(0,0) = one;
      b(0) = one;
      c(0) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(1);
    }
};


template<class Scalar>
class Implicit2Stage3rdOrderRadauB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage3rdOrderRadauB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage3rdOrderRadauB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.5, pg 74\n"
                  << "c = [ 1/3     1   ]'\n"
                  << "A = [ 5/12  -1/12 ]\n"
                  << "    [ 3/4    1/4  ]\n"
                  << "b = [ 3/4    1/4  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( 5*one/(12*one) );
      A(0,1) = as<Scalar>( -one/(12*one) );
      A(1,0) = as<Scalar>( 3*one/(4*one) );
      A(1,1) = as<Scalar>( one/(4*one) );
      b(0) = as<Scalar>( 3*one/(4*one) );
      b(1) = as<Scalar>( one/(4*one) );
      c(0) = as<Scalar>( one/(3*one) );
      c(1) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(3);
    }
};


template<class Scalar>
class Implicit3Stage5thOrderRadauB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage5thOrderRadauB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage5thOrderRadauB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.6, pg 74\n"
                  << "c = [ (4-sqrt(6))/10          (4+sqrt(6))/10                 1            ]'\n"
                  << "A = [ (88-7*sqrt(6))/360      (296-169*sqrt(6))/1800  (-2+3*sqrt(6))/225  ]\n"
                  << "    [ (296+169*sqrt(6))/1800  (88+7*sqrt(6))/360      (-2-3*sqrt(6))/225  ]\n"
                  << "    [ (16-sqrt(6))/36         (16+sqrt(6))/36         1/9                 ]\n"
                  << "b = [ (16-sqrt(6))/36         (16+sqrt(6))/36         1/9                 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( (88*one-7*one*ST::squareroot(6*one))/(360*one) );
      A(0,1) = as<Scalar>( (296*one-169*one*ST::squareroot(6*one))/(1800*one) );
      A(0,2) = as<Scalar>( (-2*one+3*one*ST::squareroot(6*one))/(225*one) );
      A(1,0) = as<Scalar>( (296*one+169*one*ST::squareroot(6*one))/(1800*one) );
      A(1,1) = as<Scalar>( (88*one+7*one*ST::squareroot(6*one))/(360*one) );
      A(1,2) = as<Scalar>( (-2*one-3*one*ST::squareroot(6*one))/(225*one) );
      A(2,0) = as<Scalar>( (16*one-ST::squareroot(6*one))/(36*one) );
      A(2,1) = as<Scalar>( (16*one+ST::squareroot(6*one))/(36*one) );
      A(2,2) = as<Scalar>( one/(9*one) );
      b(0) = as<Scalar>( (16*one-ST::squareroot(6*one))/(36*one) );
      b(1) = as<Scalar>( (16*one+ST::squareroot(6*one))/(36*one) );
      b(2) = as<Scalar>( one/(9*one) );
      c(0) = as<Scalar>( (4*one-ST::squareroot(6*one))/(10*one) );
      c(1) = as<Scalar>( (4*one+ST::squareroot(6*one))/(10*one) );
      c(2) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(5);
    }
};


template<class Scalar>
class Implicit2Stage2ndOrderLobattoA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage2ndOrderLobattoA_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage2ndOrderLobattoA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.7, pg 75\n"
                  << "c = [  0    1   ]'\n"
                  << "A = [  0    0   ]\n"
                  << "    [ 1/2  1/2  ]\n"
                  << "b = [ 1/2  1/2  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = zero;
      A(0,1) = zero;
      A(1,0) = as<Scalar>( one/(2*one) );
      A(1,1) = as<Scalar>( one/(2*one) );
      b(0) = as<Scalar>( one/(2*one) );
      b(1) = as<Scalar>( one/(2*one) );
      c(0) = zero;
      c(1) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(2);
    }
};


template<class Scalar>
class Implicit3Stage4thOrderLobattoA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage4thOrderLobattoA_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage4thOrderLobattoA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.7, pg 75\n"
                  << "c = [  0    1/2    1  ]'\n"
                  << "A = [  0     0     0   ]\n"
                  << "    [ 5/24  1/3  -1/24  ]\n"
                  << "    [ 1/6   2/3   1/6   ]\n"
                  << "b = [ 1/6   2/3   1/6   ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;
      A(1,0) = as<Scalar>( (5*one)/(24*one) );
      A(1,1) = as<Scalar>( (one)/(3*one) );
      A(1,2) = as<Scalar>( (-one)/(24*one) );
      A(2,0) = as<Scalar>( (one)/(6*one) );
      A(2,1) = as<Scalar>( (2*one)/(3*one) );
      A(2,2) = as<Scalar>( (1*one)/(6*one) );
      b(0) = as<Scalar>( (one)/(6*one) );
      b(1) = as<Scalar>( (2*one)/(3*one) );
      b(2) = as<Scalar>( (1*one)/(6*one) );
      c(0) = zero;
      c(1) = as<Scalar>( one/(2*one) );
      c(2) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Implicit4Stage6thOrderLobattoA_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit4Stage6thOrderLobattoA_RKBT()
    { 

      using Teuchos::as;
      typedef Teuchos::ScalarTraits<Scalar> ST;

      std::ostringstream description;
      description << Implicit4Stage6thOrderLobattoA_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.8, pg 75\n"
                  << "c = [ 0               (5-sqrt(5))/10       (5+sqrt(5))/10       1                 ]'\n"
                  << "A = [ 0               0                    0                    0                 ]\n"
                  << "    [ (11+sqrt(5)/120 (25-sqrt(5))/120     (25-13*sqrt(5))/120  (-1+sqrt(5))/120  ]\n"
                  << "    [ (11-sqrt(5)/120 (25+13*sqrt(5))/120  (25+sqrt(5))/120     (-1-sqrt(5))/120  ]\n"
                  << "    [ 1/12            5/12                 5/12                 1/12              ]\n"
                  << "b = [ 1/12            5/12                 5/12                 1/12              ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = zero;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;
      A(1,0) = as<Scalar>( (11*one+ST::squareroot(5*one))/(120*one) );
      A(1,1) = as<Scalar>( (25*one-ST::squareroot(5*one))/(120*one) );
      A(1,2) = as<Scalar>( (25*one-13*one*ST::squareroot(5*one))/(120*one) );
      A(1,3) = as<Scalar>( (-one+ST::squareroot(5*one))/(120*one) );
      A(2,0) = as<Scalar>( (11*one-ST::squareroot(5*one))/(120*one) );
      A(2,1) = as<Scalar>( (25*one+13*one*ST::squareroot(5*one))/(120*one) );
      A(2,2) = as<Scalar>( (25*one+ST::squareroot(5*one))/(120*one) );
      A(2,3) = as<Scalar>( (-one-ST::squareroot(5*one))/(120*one) );
      A(3,0) = as<Scalar>( one/(12*one) );
      A(3,1) = as<Scalar>( 5*one/(12*one) );
      A(3,2) = as<Scalar>( 5*one/(12*one) );
      A(3,3) = as<Scalar>( one/(12*one) );
      b(0) = as<Scalar>( one/(12*one) );
      b(1) = as<Scalar>( 5*one/(12*one) );
      b(2) = as<Scalar>( 5*one/(12*one) );
      b(3) = as<Scalar>( one/(12*one) );
      c(0) = zero;
      c(1) = as<Scalar>( (5*one-ST::squareroot(5))/(10*one) );
      c(2) = as<Scalar>( (5*one+ST::squareroot(5))/(10*one) );
      c(3) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(6);
    }
};


template<class Scalar>
class Implicit2Stage2ndOrderLobattoB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage2ndOrderLobattoB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage2ndOrderLobattoB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.9, pg 76\n"
                  << "c = [  0    1   ]'\n"
                  << "A = [ 1/2   0   ]\n"
                  << "    [ 1/2   0   ]\n"
                  << "b = [ 1/2  1/2  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(2*one) );
      A(0,1) = zero;
      A(1,0) = as<Scalar>( one/(2*one) );
      A(1,1) = zero;
      b(0) = as<Scalar>( one/(2*one) );
      b(1) = as<Scalar>( one/(2*one) );
      c(0) = zero;
      c(1) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(2);
    }
};


template<class Scalar>
class Implicit3Stage4thOrderLobattoB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage4thOrderLobattoB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage4thOrderLobattoB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.9, pg 76\n"
                  << "c = [  0    1/2    1   ]'\n"
                  << "A = [ 1/6  -1/6    0   ]\n"
                  << "    [ 1/6   1/3    0   ]\n"
                  << "    [ 1/6   5/6    0   ]\n"
                  << "b = [ 1/6   2/3   1/6  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(6*one) );
      A(0,1) = as<Scalar>( -one/(6*one) );
      A(0,2) = zero;
      A(1,0) = as<Scalar>( one/(6*one) );
      A(1,1) = as<Scalar>( one/(3*one) );
      A(1,2) = zero;
      A(2,0) = as<Scalar>( one/(6*one) );
      A(2,1) = as<Scalar>( 5*one/(6*one) );
      A(2,2) = zero;
      b(0) = as<Scalar>( one/(6*one) );
      b(1) = as<Scalar>( 2*one/(3*one) );
      b(2) = as<Scalar>( one/(6*one) );
      c(0) = zero;
      c(1) = as<Scalar>( one/(2*one) );
      c(2) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Implicit4Stage6thOrderLobattoB_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit4Stage6thOrderLobattoB_RKBT()
    { 

      std::ostringstream description;
      description << Implicit4Stage6thOrderLobattoB_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.10, pg 76\n"
                  << "c = [ 0     (5-sqrt(5))/10       (5+sqrt(5))/10       1     ]'\n"
                  << "A = [ 1/12  (-1-sqrt(5))/24      (-1+sqrt(5))/24      0     ]\n"
                  << "    [ 1/12  (25+sqrt(5))/120     (25-13*sqrt(5))/120  0     ]\n"
                  << "    [ 1/12  (25+13*sqrt(5))/120  (25-sqrt(5))/120     0     ]\n"
                  << "    [ 1/12  (11-sqrt(5))/24      (11+sqrt(5))/24      0     ]\n"
                  << "b = [ 1/12  5/12                 5/12                 1/12  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(12*one) );
      A(0,1) = as<Scalar>( (-one-ST::squareroot(5))/(24*one) );
      A(0,2) = as<Scalar>( (-one+ST::squareroot(5))/(24*one) );
      A(0,3) = zero;
      A(1,0) = as<Scalar>( one/(12*one) );
      A(1,1) = as<Scalar>( (25*one+ST::squareroot(5))/(120*one) );
      A(1,2) = as<Scalar>( (25*one-13*one*ST::squareroot(5))/(120*one) );
      A(1,3) = zero;
      A(2,0) = as<Scalar>( one/(12*one) );
      A(2,1) = as<Scalar>( (25*one+13*one*ST::squareroot(5))/(120*one) );
      A(2,2) = as<Scalar>( (25*one-ST::squareroot(5))/(120*one) );
      A(2,3) = zero;
      A(3,0) = as<Scalar>( one/(12*one) );
      A(3,1) = as<Scalar>( (11*one-ST::squareroot(5*one))/(24*one) );
      A(3,2) = as<Scalar>( (11*one+ST::squareroot(5*one))/(24*one) );
      A(3,3) = zero;
      b(0) = as<Scalar>( one/(12*one) );
      b(1) = as<Scalar>( 5*one/(12*one) );
      b(2) = as<Scalar>( 5*one/(12*one) );
      b(3) = as<Scalar>( one/(12*one) );
      c(0) = zero;
      c(1) = as<Scalar>( (5*one-ST::squareroot(5*one))/(10*one) );
      c(2) = as<Scalar>( (5*one+ST::squareroot(5*one))/(10*one) );
      c(3) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(6);
    }
};


template<class Scalar>
class Implicit2Stage2ndOrderLobattoC_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit2Stage2ndOrderLobattoC_RKBT()
    { 

      std::ostringstream description;
      description << Implicit2Stage2ndOrderLobattoC_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.11, pg 76\n"
                  << "c = [  0    1   ]'\n"
                  << "A = [ 1/2 -1/2  ]\n"
                  << "    [ 1/2  1/2  ]\n"
                  << "b = [ 1/2  1/2  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 2;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(2*one) );
      A(0,1) = as<Scalar>( -one/(2*one) );
      A(1,0) = as<Scalar>( one/(2*one) );
      A(1,1) = as<Scalar>( one/(2*one) );
      b(0) = as<Scalar>( one/(2*one) );
      b(1) = as<Scalar>( one/(2*one) );
      c(0) = zero;
      c(1) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(2);
    }
};


template<class Scalar>
class Implicit3Stage4thOrderLobattoC_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit3Stage4thOrderLobattoC_RKBT()
    { 

      std::ostringstream description;
      description << Implicit3Stage4thOrderLobattoC_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.11, pg 76\n"
                  << "c = [  0    1/2    1   ]'\n"
                  << "A = [ 1/6  -1/3   1/6  ]\n"
                  << "    [ 1/6   5/12 -1/12 ]\n"
                  << "    [ 1/6   2/3   1/6  ]\n"
                  << "b = [ 1/6   2/3   1/6  ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(6*one) );
      A(0,1) = as<Scalar>( -one/(3*one) );
      A(0,2) = as<Scalar>( one/(6*one) );
      A(1,0) = as<Scalar>( one/(6*one) );
      A(1,1) = as<Scalar>( 5*one/(12*one) );
      A(1,2) = as<Scalar>( -one/(12*one) );
      A(2,0) = as<Scalar>( one/(6*one) );
      A(2,1) = as<Scalar>( 2*one/(3*one) );
      A(2,2) = as<Scalar>( one/(6*one) );
      b(0) = as<Scalar>( one/(6*one) );
      b(1) = as<Scalar>( 2*one/(3*one) );
      b(2) = as<Scalar>( one/(6*one) );
      c(0) = zero;
      c(1) = as<Scalar>( one/(2*one) );
      c(2) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class Implicit4Stage6thOrderLobattoC_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    Implicit4Stage6thOrderLobattoC_RKBT()
    { 

      std::ostringstream description;
      description << Implicit4Stage6thOrderLobattoC_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "Table 5.12, pg 76\n"
                  << "c = [ 0     (5-sqrt(5))/10       (5+sqrt(5))/10       1          ]'\n"
                  << "A = [ 1/12  -sqrt(5)/12          sqrt(5)/12          -1/12       ]\n"
                  << "    [ 1/12  1/4                  (10-7*sqrt(5))/60   sqrt(5)/60  ]\n"
                  << "    [ 1/12  (10+7*sqrt(5))/60    1/4                 -sqrt(5)/60 ]\n"
                  << "    [ 1/12  5/12                 5/12                 1/12       ]\n"
                  << "b = [ 1/12  5/12                 5/12                 1/12       ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 4;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      A(0,0) = as<Scalar>( one/(12*one) );
      A(0,1) = as<Scalar>( -ST::squareroot(5*one)/(12*one) );
      A(0,2) = as<Scalar>( ST::squareroot(5*one)/(12*one) );
      A(0,3) = as<Scalar>( -one/(12*one) );
      A(1,0) = as<Scalar>( one/(12*one) );
      A(1,1) = as<Scalar>( one/(4*one) );
      A(1,2) = as<Scalar>( (10*one-7*one*ST::squareroot(5*one))/(60*one) );
      A(1,3) = as<Scalar>( ST::squareroot(5*one)/(60*one) );
      A(2,0) = as<Scalar>( one/(12*one) );
      A(2,1) = as<Scalar>( (10*one+7*one*ST::squareroot(5*one))/(60*one) );
      A(2,2) = as<Scalar>( one/(4*one) );
      A(2,3) = as<Scalar>( -ST::squareroot(5*one)/(60*one) );
      A(3,0) = as<Scalar>( one/(12*one) );
      A(3,1) = as<Scalar>( 5*one/(12*one) );
      A(3,2) = as<Scalar>( 5*one/(12*one) );
      A(3,3) = as<Scalar>( one/(12*one) );
      b(0) = as<Scalar>( one/(12*one) );
      b(1) = as<Scalar>( 5*one/(12*one) );
      b(2) = as<Scalar>( 5*one/(12*one) );
      b(3) = as<Scalar>( one/(12*one) );
      c(0) = zero;
      c(1) = as<Scalar>( (5*one-ST::squareroot(5*one))/(10*one) );
      c(2) = as<Scalar>( (5*one+ST::squareroot(5*one))/(10*one) );
      c(3) = one;
      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(6);
    }
};



template<class Scalar>
class SDIRK5Stage5thOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    SDIRK5Stage5thOrder_RKBT()
    { 

      std::ostringstream description;
      description << SDIRK5Stage5thOrder_name() << "\n"
        << "A-stable\n"
        << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
        << "E. Hairer and G. Wanner\n"
        << "pg101 \n"
        << "c = [ (6-sqrt(6))/10                (6+9*sqrt(6))/35             1                        (4-sqrt(6))/10         (4+sqrt(6))/10  ]'\n"
        << "A = [ (6-sqrt(6))/10                                                                                                             ]\n"
        << "    [ (-6+5*sqrt(6))/14             (6-sqrt(6))/10                                                                               ]\n"
        << "    [ (888+607*sqrt(6))/2850        (126-161*sqrt(6))/1425       (6-sqrt(6))/10                                                  ]\n"
        << "    [ (3153-3082*sqrt(6))/14250     (3213+1148*sqrt(6))/28500    (-267+88*sqrt(6))/500    (6-sqrt(6))/10                         ]\n"
        << "    [ (-32583+14638*sqrt(6))/71250  (-17199+364*sqrt(6))/142500  (1329-544*sqrt(6))/2500  (-96+131*sqrt(6))/625  (6-sqrt(6))/10  ]\n"
        << "b = [ 0                             0                            1/9                      (16-sqrt(6))/36        (16+sqrt(6))/36 ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 5;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      Scalar sqrt6 = ST::squareroot(as<Scalar>(6*one));
      Scalar gamma = as<Scalar>( (6*one - sqrt6) / (10*one) ); // diagonal
      A(0,0) = gamma;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;
      A(0,4) = zero;

      A(1,0) = as<Scalar>( (-6*one+5*one*sqrt6)/(14*one) );
      A(1,1) = gamma;
      A(1,2) = zero;
      A(1,3) = zero;
      A(1,4) = zero;

      A(2,0) = as<Scalar>( (888*one+607*one*sqrt6)/(2850*one) );
      A(2,1) = as<Scalar>( (126*one-161*one*sqrt6)/(1425*one) );
      A(2,2) = gamma;
      A(2,3) = zero;
      A(2,4) = zero;

      A(3,0) = as<Scalar>( (3153*one-3082*one*sqrt6)/(14250*one) );
      A(3,1) = as<Scalar>( (3213*one+1148*one*sqrt6)/(28500*one) );
      A(3,2) = as<Scalar>( (-267*one+88*one*sqrt6)/(500*one) );
      A(3,3) = gamma;
      A(3,4) = zero;

      A(4,0) = as<Scalar>( (-32583*one+14638*one*sqrt6)/(71250*one) );
      A(4,1) = as<Scalar>( (-17199*one+364*one*sqrt6)/(142500*one) );
      A(4,2) = as<Scalar>( (1329*one-544*one*sqrt6)/(2500*one) );
      A(4,3) = as<Scalar>( (-96*one+131*sqrt6)/(625*one) );
      A(4,4) = gamma;

      b(0) = zero;
      b(1) = zero;
      b(2) = as<Scalar>( one/(9*one) );
      b(3) = as<Scalar>( (16*one-sqrt6)/(36*one) );
      b(4) = as<Scalar>( (16*one+sqrt6)/(36*one) );

      c(0) = gamma;
      c(1) = as<Scalar>( (6*one+9*one*sqrt6)/(35*one) );
      c(2) = one;
      c(3) = as<Scalar>( (4*one-sqrt6)/(10*one) );
      c(4) = as<Scalar>( (4*one+sqrt6)/(10*one) );

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(5);
    }
};


template<class Scalar>
class SDIRK5Stage4thOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    SDIRK5Stage4thOrder_RKBT()
    { 

      std::ostringstream description;
      description << SDIRK5Stage4thOrder_name() << "\n"
        << "L-stable\n"
        << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
        << "E. Hairer and G. Wanner\n"
        << "pg100 \n"
        << "c  = [ 1/4       3/4        11/20   1/2     1   ]'\n"
        << "A  = [ 1/4                                      ]\n"
        << "     [ 1/2       1/4                            ]\n"
        << "     [ 17/50     -1/25      1/4                 ]\n"
        << "     [ 371/1360  -137/2720  15/544  1/4         ]\n"
        << "     [ 25/24     -49/48     125/16  -85/12  1/4 ]\n"
        << "b  = [ 25/24     -49/48     125/16  -85/12  1/4 ]'\n"
        << "b' = [ 59/48     -17/96     225/32  -85/12  0   ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 5;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      Scalar onequarter = as<Scalar>( one/(4*one) );
      A(0,0) = onequarter;
      A(0,1) = zero;
      A(0,2) = zero;
      A(0,3) = zero;
      A(0,4) = zero;

      A(1,0) = as<Scalar>( one / (2*one) );
      A(1,1) = onequarter;
      A(1,2) = zero;
      A(1,3) = zero;
      A(1,4) = zero;

      A(2,0) = as<Scalar>( 17*one/(50*one) );
      A(2,1) = as<Scalar>( -one/(25*one) );
      A(2,2) = onequarter;
      A(2,3) = zero;
      A(2,4) = zero;

      A(3,0) = as<Scalar>( 371*one/(1360*one) );
      A(3,1) = as<Scalar>( -137*one/(2720*one) );
      A(3,2) = as<Scalar>( 15*one/(544*one) );
      A(3,3) = onequarter;
      A(3,4) = zero;

      A(4,0) = as<Scalar>( 25*one/(24*one) );
      A(4,1) = as<Scalar>( -49*one/(48*one) );
      A(4,2) = as<Scalar>( 125*one/(16*one) );
      A(4,3) = as<Scalar>( -85*one/(12*one) );
      A(4,4) = onequarter;

      b(0) = as<Scalar>( 25*one/(24*one) );
      b(1) = as<Scalar>( -49*one/(48*one) );
      b(2) = as<Scalar>( 125*one/(16*one) );
      b(3) = as<Scalar>( -85*one/(12*one) );
      b(4) = onequarter;

      /*
      // Alternate version 
      b(0) = as<Scalar>( 59*one/(48*one) );
      b(1) = as<Scalar>( -17*one/(96*one) );
      b(2) = as<Scalar>( 225*one/(32*one) );
      b(3) = as<Scalar>( -85*one/(12*one) );
      b(4) = zero;
      */
      c(0) = onequarter;
      c(1) = as<Scalar>( 3*one/(4*one) );
      c(2) = as<Scalar>( 11*one/(20*one) );
      c(3) = as<Scalar>( one/(2*one) );
      c(4) = one;

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


template<class Scalar>
class SDIRK3Stage4thOrder_RKBT :
  virtual public RKButcherTableauDefaultBase<Scalar>
{
  public:
    SDIRK3Stage4thOrder_RKBT()
    { 

      std::ostringstream description;
      description << SDIRK3Stage4thOrder_name() << "\n"
                  << "A-stable\n"
                  << "Solving Ordinary Differential Equations II:  Stiff and Differential-Algebraic Problems, 2nd Revised Edition\n"
                  << "E. Hairer and G. Wanner\n"
                  << "pg100 \n"
                  << "gamma = (1/sqrt(3))*cos(pi/18)+1/2\n"
                  << "delta = 1/(6*(2*gamma-1)^2)\n"
                  << "c = [ gamma      1/2        1-gamma ]'\n"
                  << "A = [ gamma                         ]\n"
                  << "    [ 1/2-gamma  gamma              ]\n"
                  << "    [ 2*gamma    1-4*gamma  gamma   ]\n"
                  << "b = [ delta      1-2*delta  delta   ]'" << std::endl;
      typedef ScalarTraits<Scalar> ST;
      int numStages = 3;
      Teuchos::SerialDenseMatrix<int,Scalar> A(numStages,numStages);
      Teuchos::SerialDenseVector<int,Scalar> b(numStages);
      Teuchos::SerialDenseVector<int,Scalar> c(numStages);
      Scalar zero = ST::zero();
      Scalar one = ST::one();
      Scalar pi = as<Scalar>(4*one)*std::atan(one);
      Scalar gamma = as<Scalar>( one/ST::squareroot(3*one)*std::cos(pi/(18*one))+one/(2*one) );
      Scalar delta = as<Scalar>( one/(6*one*std::pow(2*gamma-one,2*one)) );
      A(0,0) = gamma;
      A(0,1) = zero;
      A(0,2) = zero;

      A(1,0) = as<Scalar>( one/(2*one) - gamma );
      A(1,1) = gamma;
      A(1,2) = zero;

      A(2,0) = as<Scalar>( 2*gamma );
      A(2,1) = as<Scalar>( one - 4*gamma );
      A(2,2) = gamma;

      b(0) = delta;
      b(1) = as<Scalar>( one-2*delta );
      b(2) = delta;

      c(0) = gamma;
      c(1) = as<Scalar>( one/(2*one) );
      c(2) = as<Scalar>( one - gamma );

      this->setMyDescription(description.str());
      this->setMy_A(A);
      this->setMy_b(b);
      this->setMy_c(c);
      this->setMy_order(4);
    }
};


} // namespace Rythmos


#endif // RYTHMOS_RK_BUTCHER_TABLEAU_HPP
