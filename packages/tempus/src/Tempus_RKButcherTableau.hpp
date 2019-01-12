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
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Thyra_MultiVectorStdOps.hpp"


//#include "Thyra_ProductVectorBase.hpp"

namespace Tempus {


/** \brief Runge-Kutta methods.
 *
 *  This base class specifies the Butcher tableau which defines the
 *  Runge-Kutta (RK) method.  Both explicit and implicit RK methods
 *  can be specied here, and of arbitrary number of stages and orders.
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
  virtual public Teuchos::ParameterListAcceptorDefaultBase,
  virtual public Teuchos::VerboseObject<RKButcherTableau<Scalar> >
{
  public:
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

    virtual void initialize(
      const Teuchos::SerialDenseMatrix<int,Scalar>& A,
      const Teuchos::SerialDenseVector<int,Scalar>& b,
      const Teuchos::SerialDenseVector<int,Scalar>& c,
      const int order,
      const std::string& longDescription,
      bool isEmbedded = false,
      const Teuchos::SerialDenseVector<int,Scalar>&
        bstar =  Teuchos::SerialDenseVector<int,Scalar>())
    {
      this->initialize(A,b,c,order,order,order,
                       longDescription,isEmbedded,bstar);
    }

    virtual void initialize(
      const Teuchos::SerialDenseMatrix<int,Scalar>& A,
      const Teuchos::SerialDenseVector<int,Scalar>& b,
      const Teuchos::SerialDenseVector<int,Scalar>& c,
      const int order,
      const int orderMin,
      const int orderMax,
      const std::string& longDescription,
      bool isEmbedded = false,
      const Teuchos::SerialDenseVector<int,Scalar>&
        bstar =  Teuchos::SerialDenseVector<int,Scalar>())
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
      longDescription_ = longDescription;

      if (isEmbedded) {
        TEUCHOS_ASSERT_EQUALITY( bstar.length(), numStages );
        bstar_ = bstar;
        isEmbedded_ = true;
      }
    }

    /* \brief Redefined from Teuchos::ParameterListAcceptorDefaultBase */
    //@{
      virtual void setParameterList(
        Teuchos::RCP<Teuchos::ParameterList> const& pList)
      {
        Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
        if (pList == Teuchos::null) *pl = *(this->getValidParameters());
        else pl = pList;
        pl->validateParameters(*this->getValidParameters());
        TEUCHOS_TEST_FOR_EXCEPTION(
          pl->get<std::string>("Stepper Type") != this->description()
          ,std::runtime_error,
          "  Stepper Type != \""+this->description()+"\"\n"
          "  Stepper Type = " + pl->get<std::string>("Stepper Type"));
        this->setMyParamList(pl);
        this->rkbtPL_ = pl;
      }

      virtual Teuchos::RCP<const Teuchos::ParameterList>
      getValidParameters() const
      {
        Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
        pl->setName("Default Stepper - " + this->description());
        pl->set<std::string>("Description", this->getDescription());
        pl->set<std::string>("Stepper Type", this->description());
        pl->set<bool>("Use Embedded", false);
        return pl;
      }
    //@}

    /* \brief Redefined from Teuchos::Describable */
    //@{
      virtual std::string description() const = 0;

      virtual void describe( Teuchos::FancyOStream &out,
                             const Teuchos::EVerbosityLevel verbLevel) const
      {
        if (verbLevel != Teuchos::VERB_NONE) {
          out << this->description() << std::endl;
          out << this->getDescription() << std::endl;
          out << "number of Stages = " << this->numStages() << std::endl;
          out << "A = " << this->A() << std::endl;
          out << "b = " << this->b() << std::endl;
          out << "c = " << this->c() << std::endl;
          out << "bstar = " << this->bstar() << std::endl;
          out << "order    = " << this->order()    << std::endl;
          out << "orderMin = " << this->orderMin() << std::endl;
          out << "orderMax = " << this->orderMax() << std::endl;
          out << "isImplicit = " << this->isImplicit() << std::endl;
          out << "isDIRK     = " << this->isDIRK()     << std::endl;
          out << "isEmbedded = " << this->isEmbedded() << std::endl;
        }
      }
    //@}

  protected:
    virtual void setDescription(std::string longD) { longDescription_ = longD; }
    const std::string& getDescription() const { return longDescription_; }

    void set_A(const Teuchos::SerialDenseMatrix<int,Scalar>& A) { A_ = A; }
    void set_b(const Teuchos::SerialDenseVector<int,Scalar>& b) { b_ = b; }
    void set_c(const Teuchos::SerialDenseVector<int,Scalar>& c) { c_ = c; }
    void set_order(const int& order) { order_ = order; }
    void set_orderMin(const int& order) { orderMin_ = order; }
    void set_orderMax(const int& order) { orderMax_ = order; }
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

  private:
    Teuchos::SerialDenseMatrix<int,Scalar> A_;
    Teuchos::SerialDenseVector<int,Scalar> b_;
    Teuchos::SerialDenseVector<int,Scalar> c_;
    int order_;
    int orderMin_;
    int orderMax_;
    bool isImplicit_;
    bool isDIRK_;
    std::string longDescription_;

    bool isEmbedded_ = false;
    Teuchos::SerialDenseVector<int,Scalar> bstar_;

  protected:
    Teuchos::RCP<Teuchos::ParameterList> rkbtPL_;
};

// ----------------------------------------------------------------------------
// Nonmember constructor
template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> > rKButcherTableau()
{
  return(rcp(new RKButcherTableau<Scalar>()));
}

// Nonmember constructor
template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> > rKButcherTableau(
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    int order,
    const std::string& description = ""
    )
{
  Teuchos::RCP<RKButcherTableau<Scalar> > rkbt =
    rcp(new RKButcherTableau<Scalar>());
  rkbt->initialize(A,b,c,order,order,order,description);
  return(rkbt);
}


// Nonmember constructor
template<class Scalar>
Teuchos::RCP<RKButcherTableau<Scalar> > rKButcherTableau(
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    int order,
    int orderMin,
    int orderMax,
    const std::string& description = ""
    )
{
  Teuchos::RCP<RKButcherTableau<Scalar> > rkbt =
    rcp(new RKButcherTableau<Scalar>());
  rkbt->initialize(A,b,c,order,orderMin,orderMax,description);
  return(rkbt);
}

template<class Scalar>
class General_RKButcherTableau :
  virtual public RKButcherTableau<Scalar>
{
  protected:
  void parseGeneralPL(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    using Teuchos::as;

    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., bstar).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    Teuchos::RCP<Teuchos::ParameterList> tableauPL = sublist(pl,"Tableau",true);
    std::size_t numStages = 0;
    int order = tableauPL->get<int>("order");
    Teuchos::SerialDenseMatrix<int,Scalar> A;
    Teuchos::SerialDenseVector<int,Scalar> b;
    Teuchos::SerialDenseVector<int,Scalar> c;
    Teuchos::SerialDenseVector<int,Scalar> bstar;

    // read in the A matrix
    {
      std::vector<std::string> A_row_tokens;
      Tempus::StringTokenizer(A_row_tokens, tableauPL->get<std::string>("A"),
                              ";",true);

      // this is the only place where numStages is set
      numStages = A_row_tokens.size();

      // allocate the matrix
      A.shape(as<int>(numStages),as<int>(numStages));

      // fill the rows
      for(std::size_t row=0;row<numStages;row++) {
        // parse the row (tokenize on space)
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(tokens,A_row_tokens[row]," ",true);

        std::vector<double> values;
        Tempus::TokensToDoubles(values,tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
          "Error parsing A matrix, wrong number of stages in row "
          << row << "\n" + this->description());

        for(std::size_t col=0;col<numStages;col++)
          A(row,col) = values[col];
      }
    }

    // size b and c vectors
    b.size(as<int>(numStages));
    c.size(as<int>(numStages));

    // read in the b vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,tableauPL->get<std::string>("b")," ",true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
        "Error parsing b vector, wrong number of stages.\n"
        + this->description());

      for(std::size_t i=0;i<numStages;i++)
        b(i) = values[i];
    }

    // read in the c vector
    {
      std::vector<std::string> tokens;
      Tempus::StringTokenizer(tokens,tableauPL->get<std::string>("c")," ",true);
      std::vector<double> values;
      Tempus::TokensToDoubles(values,tokens);

      TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
        "Error parsing c vector, wrong number of stages.\n"
        + this->description());

      for(std::size_t i=0;i<numStages;i++)
        c(i) = values[i];
    }

    if (tableauPL->isParameter("bstar") and
        tableauPL->get<std::string>("bstar") != "1.0") {
      bstar.size(as<int>(numStages));
      // read in the bstar vector
      {
        std::vector<std::string> tokens;
        Tempus::StringTokenizer(
          tokens, tableauPL->get<std::string>("bstar"), " ", true);
        std::vector<double> values;
        Tempus::TokensToDoubles(values,tokens);

        TEUCHOS_TEST_FOR_EXCEPTION(values.size()!=numStages,std::runtime_error,
          "Error parsing bstar vector, wrong number of stages.\n"
          + this->description());

        for(std::size_t i=0;i<numStages;i++)
          bstar(i) = values[i];
      }
      this->initialize(A,b,c,order,this->getDescription(), true, bstar);
    } else {
      this->initialize(A,b,c,order,this->getDescription());
    }

    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

};

// ----------------------------------------------------------------------------
/** \brief General Explicit Runge-Kutta Butcher Tableau
 *
 *  The format of the Butcher Tableau parameter list is
    \verbatim
      <Parameter name="A" type="string" value="# # # ;
                                               # # # ;
                                               # # #">
      <Parameter name="b" type="string" value="# # #">
      <Parameter name="c" type="string" value="# # #">
    \endverbatim
 *  Note the number of stages is implicit in the number of entries.
 *  The number of stages must be consistent.
 *
 *  Default tableau is RK4 (order=4):
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                        1/2 &  0  & 1/2 &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 1/3 & 1/3 & 1/6 \end{array}
 *  \f]
 */
template<class Scalar>
class GeneralExplicit_RKBT :
  virtual public General_RKButcherTableau<Scalar>
{
  public:
  GeneralExplicit_RKBT()
  {
    std::stringstream Description;
    Description << this->description() << "\n"
      << "The format of the Butcher Tableau parameter list is\n"
      << "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
      << "                                           # # # ;\n"
      << "                                           # # #\"/>\n"
      << "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
      << "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
      << "Note the number of stages is implicit in the number of entries.\n"
      << "The number of stages must be consistent.\n"
      << "\n"
      << "Default tableau is RK4 (order=4):\n"
      << "c = [  0  1/2 1/2  1  ]'\n"
      << "A = [  0              ]\n"
      << "    [ 1/2  0          ]\n"
      << "    [  0  1/2  0      ]\n"
      << "    [  0   0   1   0  ]\n"
      << "b = [ 1/6 1/3 1/3 1/6 ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "General ERK"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    this->parseGeneralPL(pList);
    TEUCHOS_TEST_FOR_EXCEPTION(this->isImplicit() == true, std::logic_error,
      "Error - General ERK received an implicit Butcher Tableau!\n");
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);

    // Tableau ParameterList
    Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
    tableauPL->set<std::string>("A",
     "0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.0; 0.0 0.5 0.0 0.0; 0.0 0.0 1.0 0.0");
    tableauPL->set<std::string>("b",
     "0.166666666666667 0.333333333333333 0.333333333333333 0.166666666666667");
    //tableauPL->set<std::string>("bstar", "1.0");
    tableauPL->set<std::string>("c", "0.0 0.5 0.5 1.0");
    tableauPL->set<int>("order", 4);
    pl->set("Tableau", *tableauPL);

    return pl;
  }
};


// ----------------------------------------------------------------------------
/** \brief Backward Euler Runge-Kutta Butcher Tableau
 *
 *  The tableau for Backward Euler (order=1) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} 1 & 1 \\ \hline
 *                       & 1 \end{array}
 *  \f]
 */
template<class Scalar>
class BackwardEuler_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  BackwardEuler_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "c = [ 1 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "RK Backward Euler"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    A(0,0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    b(0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    c(0) = ST::one();
    int order = 1;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};


// ----------------------------------------------------------------------------
/** \brief Forward Euler Runge-Kutta Butcher Tableau
 *
 *  The tableau for Forward Euler (order=1) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} 0 & 0 \\ \hline
 *                       & 1 \end{array}
 *  \f]
 */
template<class Scalar>
class ForwardEuler_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  ForwardEuler_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "c = [ 0 ]'\n"
                << "A = [ 0 ]\n"
                << "b = [ 1 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    b(0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    int order = 1;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const { return "RK Forward Euler"; }
};


// ----------------------------------------------------------------------------
/** \brief Runge-Kutta 4th order Butcher Tableau
 *
 *  The tableau for RK4 (order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                        1/2 &  0  & 1/2 &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 1/3 & 1/3 & 1/6 \end{array}
 *  \f]
 */
template<class Scalar>
class Explicit4Stage4thOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit4Stage4thOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "\"The\" Runge-Kutta Method (explicit):\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.2, pg 138\n"
                << "c = [  0  1/2 1/2  1  ]'\n"
                << "A = [  0              ] \n"
                << "    [ 1/2  0          ]\n"
                << "    [  0  1/2  0      ]\n"
                << "    [  0   0   1   0  ]\n"
                << "b = [ 1/6 1/3 1/3 1/6 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar onethird = one/(3*one);

    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

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

    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const { return "RK Explicit 4 Stage"; }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RK Bogacki-Shampine Butcher Tableau
 *
 *  The tableau (order=3(2)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\ \hline
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  & 0    &     &     & \\
 *                        1/3 & 1/2  & 0   &     & \\
 *                        2/3 & 0    & 3/4 & 0   & \\
 *                         1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                            & 2/9  & 1/3 & 4/9 & 0 \\
 *                            & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *  \f]
 *  Reference:  P. Bogacki and L.F. Shampine.
 *              A 3(2) pair of Runge–Kutta formulas.
 *              Applied Mathematics Letters, 2(4):321 – 325, 1989.
 *
 */
template<class Scalar>
class ExplicitBogackiShampine32_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  ExplicitBogackiShampine32_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "P. Bogacki and L.F. Shampine.\n"
                << "A 3(2) pair of Runge–Kutta formulas.\n"
                << "Applied Mathematics Letters, 2(4):321 – 325, 1989.\n"
                << "c =     [ 0     1/3  2/3   1  ]'\n"
                << "A =     [ 0                   ]\n"
                << "        [ 1/2    0            ]\n"
                << "        [  0    3/4   0       ]\n"
                << "        [ 2/9   1/3  4/9   0  ]\n"
                << "b     = [ 2/9   1/3  4/9   0  ]\n"
                << "bstar = [ 7/24  1/4  1/3  1/8 ]\n" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) = zero;
    A(0,1) = zero;
    A(0,2) = zero;
    A(0,3) = zero;

    A(1,0) = as<Scalar>(one/(2*one));
    A(1,1) = zero;
    A(1,2) = zero;
    A(1,3) = zero;

    A(2,0) = zero;
    A(2,1) = as<Scalar>(3*one/(4*one));
    A(2,2) = zero;
    A(2,3) = zero;

    A(3,0) = as<Scalar>(2*one/(9*one));
    A(3,1) = as<Scalar>(1*one/(3*one));
    A(3,2) = as<Scalar>(4*one/(9*one));
    A(3,3) = zero;

    // Fill b:
    b(0) = A(3,0);
    b(1) = A(3,1);
    b(2) = A(3,2);
    b(3) = A(3,3);

    // Fill c:
    c(0) = zero;
    c(1) = as<Scalar>(1*one/(3*one));
    c(2) = as<Scalar>(2*one/(3*one));
    c(3) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(7.0*one/(24*one));
    bstar(1) = as<Scalar>(1*one/(4*one));
    bstar(2) = as<Scalar>(1*one/(3*one));
    bstar(3) = as<Scalar>(1*one/(8*one));
    int order = 3;

    this->initialize(A,b,c,order,Description.str(),true,bstar);
  }
  virtual std::string description() const { return "Bogacki-Shampine 3(2) Pair"; }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RK Merson Butcher Tableau
 *
 *  The tableau (order=4(5)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\ \hline
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}  0 & 0    &     &      &     & \\
 *                        1/3 & 1/3  & 0   &      &     & \\
 *                        1/3 & 1/6  & 1/6 & 0    &     & \\
 *                        1/2 & 1/8  & 0   & 3/8  &     & \\
 *                         1  & 1/2  & 0   & -3/2 & 2   & \\ \hline
 *                            & 1/6  & 0   & 0    & 2/3 & 1/6 \\
 *                            & 1/10 & 0   & 3/10 & 2/5 & 1/5 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 4.1, pg 167.
 *
 */
template<class Scalar>
class ExplicitMerson45_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  ExplicitMerson45_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 4.1, pg 167\n"
                << "c =     [  0    1/3  1/3  1/2   1  ]'\n"
                << "A =     [  0                       ]\n"
                << "        [ 1/3    0                 ]\n"
                << "        [ 1/6   1/6   0            ]\n"
                << "        [ 1/8    0   3/8   0       ]\n"
                << "        [ 1/2    0  -3/2   2    0  ]\n"
                << "b     = [ 1/6    0    0   2/3  1/6 ]\n"
                << "bstar = [ 1/10   0  3/10  2/5  1/5 ]\n" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages, true);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages, true);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(1,0) = as<Scalar>(one/(3*one));;

    A(2,0) = as<Scalar>(one/(6*one));;
    A(2,1) = as<Scalar>(one/(6*one));;

    A(3,0) = as<Scalar>(one/(8*one));;
    A(3,2) = as<Scalar>(3*one/(8*one));;

    A(4,0) = as<Scalar>(one/(2*one));;
    A(4,2) = as<Scalar>(-3*one/(2*one));;
    A(4,3) = 2*one;

    // Fill b:
    b(0) = as<Scalar>(one/(6*one));
    b(3) = as<Scalar>(2*one/(3*one));
    b(4) = as<Scalar>(one/(6*one));

    // Fill c:
    c(0) = zero;
    c(1) = as<Scalar>(1*one/(3*one));
    c(2) = as<Scalar>(1*one/(3*one));
    c(3) = as<Scalar>(1*one/(2*one));
    c(4) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(1*one/(10*one));
    bstar(2) = as<Scalar>(3*one/(10*one));
    bstar(3) = as<Scalar>(2*one/(5*one));
    bstar(4) = as<Scalar>(1*one/(5*one));
    int order = 4;

    this->initialize(A,b,c,order,Description.str(),true,bstar);
  }
  virtual std::string description() const { return "Merson 4(5) Pair"; }
};

// ----------------------------------------------------------------------------
/** \brief Explicit RK 3/8th Rule Butcher Tableau
 *
 *  The tableau (order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/3 & 1/3 &  0  &     &    \\
 *                        2/3 &-1/3 &  1  &  0  &    \\
 *                         1  &  1  & -1  &  1  &  0 \\ \hline
 *                            & 1/8 & 3/8 & 3/8 & 1/8 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.2, pg 138.
 */
template<class Scalar>
class Explicit3_8Rule_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit3_8Rule_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.2, pg 138\n"
                << "c = [  0  1/3 2/3  1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/3  0          ]\n"
                << "    [-1/3  1   0      ]\n"
                << "    [  1  -1   1   0  ]\n"
                << "b = [ 1/8 3/8 3/8 1/8 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar one_third    = as<Scalar>(one/(3*one));
    const Scalar two_third    = as<Scalar>(2*one/(3*one));
    const Scalar one_eighth   = as<Scalar>(one/(8*one));
    const Scalar three_eighth = as<Scalar>(3*one/(8*one));

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

    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const { return "RK Explicit 3/8 Rule"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 4 Stage 3rd order by Runge
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  &  0  &     &     &    \\
 *                        1/2 & 1/2 &  0  &     &    \\
 *                         1  &  0  &  1  &  0  &    \\
 *                         1  &  0  &  0  &  1  &  0 \\ \hline
 *                            & 1/6 & 2/3 &  0  & 1/6 \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class Explicit4Stage3rdOrderRunge_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit4Stage3rdOrderRunge_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2  1   1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/2  0          ]\n"
                << "    [  0   1   0      ]\n"
                << "    [  0   0   1   0  ]\n"
                << "b = [ 1/6 2/3  0  1/6 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar zero = ST::zero();

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

    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 4 Stage 3rd order by Runge"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 5 Stage 3rd order by Kinnmark and Gray
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}  0  &  0  &     &     &     &    \\
 *                         1/5 & 1/5 &  0  &     &     &    \\
 *                         1/5 &  0  & 1/5 &  0  &     &    \\
 *                         1/3 &  0  &  0  & 1/3 &  0  &    \\
 *                         2/3 &  0  &  0  &  0  & 2/3 &  0 \\ \hline
 *                             & 1/4 &  0  &  0  &  0  & 3/4 \end{array}
 *  \f]
 *  Reference:  Modified by P. Ullrich.  From the prim_advance_mod.F90
 *              routine in the HOMME atmosphere model code.
 */
template<class Scalar>
class Explicit5Stage3rdOrderKandG_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit5Stage3rdOrderKandG_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Kinnmark & Gray 5 stage, 3rd order scheme \n"
                << "Modified by P. Ullrich.  From the prim_advance_mod.F90 \n"
                << "routine in the HOMME atmosphere model code.\n"
                << "c = [  0  1/5  1/5  1/3  2/3  ]'\n"
                << "A = [  0                      ]\n"
                << "    [ 1/5  0                  ]\n"
                << "    [  0  1/5   0             ]\n"
                << "    [  0   0   1/3   0        ]\n"
                << "    [  0   0    0   2/3   0   ]\n"
                << "b = [ 1/4  0    0    0   3/4  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar onefifth = one/(5*one);
    const Scalar onefourth = one/(4*one);
    const Scalar onethird = one/(3*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar threefourths = 3*one/(4*one);
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) = zero;
    A(0,1) = zero;
    A(0,2) = zero;
    A(0,3) = zero;
    A(0,4) = zero;

    A(1,0) = onefifth;
    A(1,1) = zero;
    A(1,2) = zero;
    A(1,3) = zero;
    A(1,4) = zero;

    A(2,0) = zero;
    A(2,1) = onefifth;
    A(2,2) = zero;
    A(2,3) = zero;
    A(2,4) = zero;

    A(3,0) = zero;
    A(3,1) = zero;
    A(3,2) = onethird;
    A(3,3) = zero;
    A(3,4) = zero;

    A(4,0) = zero;
    A(4,1) = zero;
    A(4,2) = zero;
    A(4,3) = twothirds;
    A(4,4) = zero;

    // Fill b:
    b(0) = onefourth;
    b(1) = zero;
    b(2) = zero;
    b(3) = zero;
    b(4) = threefourths;

    // Fill c:
    c(0) = zero;
    c(1) = onefifth;
    c(2) = onefifth;
    c(3) = onethird;
    c(4) = twothirds;

    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 5 Stage 3rd order by Kinnmark and Gray"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 3 Stage 3rd order
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                       1/2 & 1/2 &  0  &     \\
 *                        1  & -1  &  2  &  0  \\ \hline
 *                           & 1/6 & 4/6 & 1/6  \end{array}
 *  \f]
 */
template<class Scalar>
class Explicit3Stage3rdOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit3Stage3rdOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "c = [  0  1/2  1  ]'\n"
                << "A = [  0          ]\n"
                << "    [ 1/2  0      ]\n"
                << "    [ -1   2   0  ]\n"
                << "b = [ 1/6 4/6 1/6 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar two = Teuchos::as<Scalar>(2*one);
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onesixth = one/(6*one);
    const Scalar foursixth = 4*one/(6*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

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

    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 3 Stage 3rd order"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 3 Stage 3rd order TVD
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                        1  &  1  &  0  &     \\
 *                       1/2 & 1/4 & 1/4 &  0  \\ \hline
 *                           & 1/6 & 1/6 & 4/6  \end{array}
 *  \f]
 *  Reference: Sigal Gottlieb and Chi-Wang Shu,
 *             'Total Variation Diminishing Runge-Kutta Schemes',
 *             Mathematics of Computation,
 *             Volume 67, Number 221, January 1998, pp. 73-85.
 *
 *  This is also written in the following set of updates.
    \verbatim
      u1 = u^n + dt L(u^n)
      u2 = 3 u^n/4 + u1/4 + dt L(u1)/4
      u^(n+1) = u^n/3 + 2 u2/2 + 2 dt L(u2)/3
    \endverbatim
 */
template<class Scalar>
class Explicit3Stage3rdOrderTVD_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit3Stage3rdOrderTVD_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                  << "Sigal Gottlieb and Chi-Wang Shu\n"
                  << "`Total Variation Diminishing Runge-Kutta Schemes'\n"
                  << "Mathematics of Computation\n"
                  << "Volume 67, Number 221, January 1998, pp. 73-85\n"
                  << "c = [  0   1  1/2 ]'\n"
                  << "A = [  0          ]\n"
                  << "    [  1   0      ]\n"
                  << "    [ 1/4 1/4  0  ]\n"
                  << "b = [ 1/6 1/6 4/6 ]'\n"
                  << "This is also written in the following set of updates.\n"
                  << "u1 = u^n + dt L(u^n)\n"
                  << "u2 = 3 u^n/4 + u1/4 + dt L(u1)/4\n"
                  << "u^(n+1) = u^n/3 + 2 u2/2 + 2 dt L(u2)/3"
                  << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onefourth = one/(4*one);
    const Scalar onesixth = one/(6*one);
    const Scalar foursixth = 4*one/(6*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) = zero;
    A(0,1) = zero;
    A(0,2) = zero;

    A(1,0) = one;
    A(1,1) = zero;
    A(1,2) = zero;

    A(2,0) = onefourth;
    A(2,1) = onefourth;
    A(2,2) = zero;

    // Fill b:
    b(0) = onesixth;
    b(1) = onesixth;
    b(2) = foursixth;

    // fill b_c_
    c(0) = zero;
    c(1) = one;
    c(2) = onehalf;

    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 3 Stage 3rd order TVD"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 3 Stage 3rd order by Heun
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  0  &  0  &     &     \\
 *                       1/3 & 1/3 &  0  &     \\
 *                       2/3 &  0  & 2/3 &  0  \\ \hline
 *                           & 1/4 &  0  & 3/4  \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class Explicit3Stage3rdOrderHeun_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit3Stage3rdOrderHeun_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/3 2/3 ]'\n"
                << "A = [  0          ] \n"
                << "    [ 1/3  0      ]\n"
                << "    [  0  2/3  0  ]\n"
                << "b = [ 1/4  0  3/4 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onethird = one/(3*one);
    const Scalar twothirds = 2*one/(3*one);
    const Scalar onefourth = one/(4*one);
    const Scalar threefourths = 3*one/(4*one);

    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

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

    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 3 Stage 3rd order by Heun"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit 2 Stage 2nd order by Runge
 *
 *  The tableau (order=2) (also known as Explicit Midpoint) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                      1/2 & 1/2 &  0  \\ \hline
 *                          &  0  &  1   \end{array}
 *  \f]
 *  Reference:  E. Hairer, S.P. Norsett, G. Wanner,
 *              "Solving Ordinary Differential Equations I:
 *              Nonstiff Problems", 2nd Revised Edition,
 *              Table 1.1, pg 135.
 */
template<class Scalar>
class Explicit2Stage2ndOrderRunge_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Explicit2Stage2ndOrderRunge_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Also known as Explicit Midpoint\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2 ]'\n"
                << "A = [  0      ]\n"
                << "    [ 1/2  0  ]\n"
                << "b = [  0   1  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

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

    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Explicit 2 Stage 2nd order by Runge"; }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit Trapezoidal
 *
 *  The tableau (order=2) (also known as Explicit Midpoint) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                       1  &  1  &  0  \\ \hline
 *                          & 1/2 & 1/2  \end{array}
 *  \f]
 */
template<class Scalar>
class ExplicitTrapezoidal_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  ExplicitTrapezoidal_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "c = [  0   1  ]'\n"
                << "A = [  0      ]\n"
                << "    [  1   0  ]\n"
                << "b = [ 1/2 1/2 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) = zero;
    A(0,1) = zero;

    A(1,0) = one;
    A(1,1) = zero;

    // Fill b:
    b(0) = onehalf;
    b(1) = onehalf;

    // fill b_c_
    c(0) = zero;
    c(1) = one;

    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const { return "RK Explicit Trapezoidal"; }
};


// ----------------------------------------------------------------------------
/** \brief General Implicit Runge-Kutta Butcher Tableau
 *
 *  The format of the Butcher Tableau parameter list is
    \verbatim
      <Parameter name="A" type="string" value="# # # ;
                                               # # # ;
                                               # # #">
      <Parameter name="b" type="string" value="# # #">
      <Parameter name="c" type="string" value="# # #">
    \endverbatim
 *  Note the number of stages is implicit in the number of entries.
 *  The number of stages must be consistent.
 *
 *  Default tableau is "SDIRK 2 Stage 2nd order":
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc} \gamma  & \gamma &        \\
 *                         1    & 1-\gamma & \gamma \\ \hline
 *                              & 1-\gamma & \gamma  \end{array}
 *  \f]
 *  where \f$\gamma = (2\pm \sqrt{2})/2\f$.  This will produce an
 *  L-stable 2nd order method.
 *
 *  Reference: U. M. Ascher and L. R. Petzold,
 *             Computer Methods for ODEs and DAEs, p. 106.
 */
template<class Scalar>
class GeneralDIRK_RKBT :
  virtual public General_RKButcherTableau<Scalar>
{
  public:
  GeneralDIRK_RKBT()
  {
    std::stringstream Description;
    Description << this->description() << "\n"
      << "The format of the Butcher Tableau parameter list is\n"
      << "  <Parameter name=\"A\" type=\"string\" value=\"# # # ;\n"
      << "                                           # # # ;\n"
      << "                                           # # #\"/>\n"
      << "  <Parameter name=\"b\" type=\"string\" value=\"# # #\"/>\n"
      << "  <Parameter name=\"c\" type=\"string\" value=\"# # #\"/>\n\n"
      << "Note the number of stages is implicit in the number of entries.\n"
      << "The number of stages must be consistent.\n"
      << "\n"
      << "Default tableau is 'SDIRK 2 Stage 2nd order':\n"
      << "  Computer Methods for ODEs and DAEs\n"
      << "  U. M. Ascher and L. R. Petzold\n"
      << "  p. 106\n"
      << "  gamma = (2+-sqrt(2))/2\n"
      << "  c = [  gamma   1     ]'\n"
      << "  A = [  gamma   0     ]\n"
      << "      [ 1-gamma  gamma ]\n"
      << "  b = [ 1-gamma  gamma ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "General DIRK"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    this->parseGeneralPL(pList);
    TEUCHOS_TEST_FOR_EXCEPTION(this->isImplicit() != true, std::logic_error,
      "Error - General DIRK did not receive a DIRK Butcher Tableau!\n");
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);

    // Tableau ParameterList
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    std::string gamma = std::to_string(Teuchos::as<Scalar>((2*one-ST::squareroot(2*one))/(2*one)));
    std::string one_gamma = std::to_string(Teuchos::as<Scalar>(one-(2*one-ST::squareroot(2*one))/(2*one)));
    Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
    tableauPL->set<std::string>("A", gamma + " 0.0; " + one_gamma + " "+gamma);
    tableauPL->set<std::string>("b", one_gamma + " " + gamma);
    tableauPL->set<std::string>("c", gamma + " 1.0");
    tableauPL->set<int>("order", 2);
    pl->set("Tableau", *tableauPL);

    return pl;
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 1 Stage 1st order
 *
 *  The tableau (order=1) (also known as Backward Euler) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} 1 & 1 \\ \hline
 *                       & 1 \end{array}
 *  \f]
 */
template<class Scalar>
class SDIRK1Stage1stOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK1Stage1stOrder_RKBT()
  {
    std::stringstream Description;
    Description << this->description() << "\n"
                << "c = [ 1 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 1 Stage 1st order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    A(0,0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    b(0) = ST::one();
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    c(0) = ST::one();
    int order = 1;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 2 Stage 2nd order
 *
 *  The tableau (order=1 or 2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc} \gamma  & \gamma &        \\
 *                         1    & 1-\gamma & \gamma \\ \hline
 *                              & 1-\gamma & \gamma  \end{array}
 *  \f]
 *  The default value is \f$\gamma = (2\pm \sqrt{2})/2\f$.
 *  This will produce an L-stable 2nd order method with the stage
 *  times within the timestep.  Other values of gamma will still
 *  produce an L-stable scheme, but will only be 1st order accurate.
 *  L-stability is guaranteed because \f$A_{sj} = b_j\f$.
 *
 *  Reference: U. M. Ascher and L. R. Petzold,
 *             Computer Methods for ODEs and DAEs, p. 106.
 */
template<class Scalar>
class SDIRK2Stage2ndOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK2Stage2ndOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Computer Methods for ODEs and DAEs\n"
                << "U. M. Ascher and L. R. Petzold\n"
                << "p. 106\n"
                << "gamma = (2+-sqrt(2))/2\n"
                << "c = [  gamma   1     ]'\n"
                << "A = [  gamma   0     ]\n"
                << "    [ 1-gamma  gamma ]\n"
                << "b = [ 1-gamma  gamma ]'" << std::endl;

    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    gamma_default_ = Teuchos::as<Scalar>((2*one-ST::squareroot(2*one))/(2*one));
    gamma_ = gamma_default_;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 2 Stage 2nd order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    gamma_ = pl->get<double>("gamma", gamma_default_);

    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    A(0,0) = gamma_;
    A(0,1) = zero;
    A(1,0) = Teuchos::as<Scalar>( one - gamma_ );
    A(1,1) = gamma_;
    b(0) = Teuchos::as<Scalar>( one - gamma_ );
    b(1) = gamma_;
    c(0) = gamma_;
    c(1) = one;

    int order = 1;
    if ( std::abs((gamma_-gamma_default_)/gamma_) < 1.0e-08 ) order = 2;

    this->initialize(A,b,c,order,1,2,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");
    pl->set<double>("gamma",gamma_default_,
      "The default value is gamma = (2-sqrt(2))/2. "
      "This will produce an L-stable 2nd order method with the stage "
      "times within the timestep.  Other values of gamma will still "
      "produce an L-stable scheme, but will only be 1st order accurate.");

    return pl;
  }

  private:
    Scalar gamma_default_;
    Scalar gamma_;
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 2 Stage 3rd order
 *
 *  The tableau (order=2 or 3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  \gamma  &  \gamma   &        \\
 *                      1-\gamma & 1-2\gamma & \gamma \\ \hline
 *                               &   1/2     &   1/2   \end{array}
 *  \f]
 *  \f[
 *  \gamma = \left\{ \begin{array}{cc}
 *                     (2\pm \sqrt{2})/2 & \mbox{then 2nd order and L-stable} \\
 *                     (3\pm \sqrt{3})/6 & \mbox{then 3rd order and A-stable}
 *                   \end{array} \right.
 *  \f]
 *  The default value is \f$\gamma = (3\pm \sqrt{3})/6\f$.
 *
 *  Reference: E. Hairer, S. P. Norsett, and G. Wanner,
 *             Solving Ordinary Differential Equations I:
 *             Nonstiff Problems, 2nd Revised Edition,
 *             Table 7.2, pg 207.
 */
template<class Scalar>
class SDIRK2Stage3rdOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK2Stage3rdOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.2, pg 207\n"
                << "gamma = (3+-sqrt(3))/6 -> 3rd order and A-stable\n"
                << "gamma = (2+-sqrt(2))/2 -> 2nd order and L-stable\n"
                << "c = [  gamma     1-gamma  ]'\n"
                << "A = [  gamma     0        ]\n"
                << "    [ 1-2*gamma  gamma    ]\n"
                << "b = [ 1/2        1/2      ]'" << std::endl;

    thirdOrderAStable_default_ = true;
    secondOrderLStable_default_ = false;
    thirdOrderAStable_ = true;
    secondOrderLStable_ = false;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    gamma_default_ =
      Teuchos::as<Scalar>( (3*one + ST::squareroot(3*one))/(6*one) );
    gamma_ = gamma_default_;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., bstar).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    thirdOrderAStable_  = pl->get<bool>("3rd Order A-stable",false);
    secondOrderLStable_ = pl->get<bool>("2nd Order L-stable",false);
    TEUCHOS_TEST_FOR_EXCEPTION(
      thirdOrderAStable_ && secondOrderLStable_, std::logic_error,
      "'3rd Order A-stable' and '2nd Order L-stable' can not both be true.");
    gamma_ = pl->get<double>("gamma", gamma_default_);

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar gammaLStable =
      as<Scalar>( (2*one + ST::squareroot(2*one))/(2*one) );
    if (thirdOrderAStable_)
      gamma_ = gamma_default_;
    else if (secondOrderLStable_)
      gamma_ = gammaLStable;
    A(0,0) = gamma_;
    A(0,1) = zero;
    A(1,0) = as<Scalar>( one - 2*gamma_ );
    A(1,1) = gamma_;
    b(0) = as<Scalar>( one/(2*one) );
    b(1) = as<Scalar>( one/(2*one) );
    c(0) = gamma_;
    c(1) = as<Scalar>( one - gamma_ );

    int order = 2;
    if ( std::abs((gamma_-gamma_default_)/gamma_) < 1.0e-08 ) {
      order = 3;
      thirdOrderAStable_  = true;
      secondOrderLStable_ = false;
    } else if ( std::abs((gamma_-gammaLStable)/gamma_) < 1.0e-08 ) {
      thirdOrderAStable_  = false;
      secondOrderLStable_ = true;
    } else {
      thirdOrderAStable_  = false;
      secondOrderLStable_ = false;
    }

    this->initialize(A,b,c,order,2,3,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");
    pl->set<bool>("3rd Order A-stable",thirdOrderAStable_default_,
      "If true, set gamma to gamma = (3+sqrt(3))/6 to obtain "
      "a 3rd order A-stable scheme. '3rd Order A-stable' and "
      "'2nd Order L-stable' can not both be true.");
    pl->set<bool>("2nd Order L-stable",secondOrderLStable_default_,
      "If true, set gamma to gamma = (2+sqrt(2))/2 to obtain "
      "a 2nd order L-stable scheme. '3rd Order A-stable' and "
      "'2nd Order L-stable' can not both be true.");
    pl->set<double>("gamma",gamma_default_,
      "If both '3rd Order A-stable' and '2nd Order L-stable' "
      "are false, gamma will be used. The default value is the "
      "'3rd Order A-stable' gamma value, (3+sqrt(3))/6.");

    return pl;
  }

  virtual std::string description() const { return "SDIRK 2 Stage 3rd order"; }

  private:
    bool thirdOrderAStable_default_;
    bool thirdOrderAStable_;
    bool secondOrderLStable_default_;
    bool secondOrderLStable_;
    Scalar gamma_default_;
    Scalar gamma_;
};


// ----------------------------------------------------------------------------
/** \brief EDIRK 2 Stage 3rd order
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                      2/3 & 1/3 & 1/3 \\ \hline
 *                          & 1/4 & 3/4  \end{array}
 *  \f]
 *  Reference: E. Hairer, S. P. Norsett, and G. Wanner,
 *             Solving Ordinary Differential Equations I:
 *             Nonstiff Problems, 2nd Revised Edition,
 *             Table 7.1, pg 205.
 */
template<class Scalar>
class EDIRK2Stage3rdOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  EDIRK2Stage3rdOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Hammer & Hollingsworth method\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.1, pg 205\n"
                << "c = [  0   2/3 ]'\n"
                << "A = [  0    0  ]\n"
                << "    [ 1/3  1/3 ]\n"
                << "b = [ 1/4  3/4 ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "EDIRK 2 Stage 3rd order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));
    this->setDescription(pl->get<std::string>("Description",""));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    A(0,0) = zero;
    A(0,1) = zero;
    A(1,0) = as<Scalar>( one/(3*one) );
    A(1,1) = as<Scalar>( one/(3*one) );
    b(0) = as<Scalar>( one/(4*one) );
    b(1) = as<Scalar>( 3*one/(4*one) );
    c(0) = zero;
    c(1) = as<Scalar>( 2*one/(3*one) );
    int order = 3;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }

};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage6thOrderKuntzmannButcher_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage6thOrderKuntzmannButcher_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "Kuntzmann & Butcher method\n"
      << "Solving Ordinary Differential Equations I:\n"
      << "Nonstiff Problems, 2nd Revised Edition\n"
      << "E. Hairer, S. P. Norsett, and G. Wanner\n"
      << "Table 7.4, pg 209\n"
      << "c = [ 1/2-sqrt(15)/10   1/2              1/2+sqrt(15)/10  ]'\n"
      << "A = [ 5/36              2/9-sqrt(15)/15  5/36-sqrt(15)/30 ]\n"
      << "    [ 5/36+sqrt(15)/24  2/9              5/36-sqrt(15)/24 ]\n"
      << "    [ 5/36+sqrt(15)/30  2/9+sqrt(15)/15  5/36             ]\n"
      << "b = [ 5/18              4/9              5/18             ]'"
      << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
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
    int order = 6;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 6th Order Kuntzmann & Butcher"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit4Stage8thOrderKuntzmannButcher_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit4Stage8thOrderKuntzmannButcher_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Kuntzmann & Butcher method\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
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
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar onehalf = as<Scalar>( one/(2*one) );
    const Scalar w1 = as<Scalar>( one/(8*one) - ST::squareroot(30*one)/(144*one) );
    const Scalar w2 = as<Scalar>( (one/(2*one))*ST::squareroot((15*one+2*one*ST::squareroot(30*one))/(35*one)) );
    const Scalar w3 = as<Scalar>( w2*(one/(6*one)+ST::squareroot(30*one)/(24*one)) );
    const Scalar w4 = as<Scalar>( w2*(one/(21*one)+5*one*ST::squareroot(30*one)/(168*one)) );
    const Scalar w5 = as<Scalar>( w2-2*w3 );
    const Scalar w1p = as<Scalar>( one/(8*one) + ST::squareroot(30*one)/(144*one) );
    const Scalar w2p = as<Scalar>( (one/(2*one))*ST::squareroot((15*one-2*one*ST::squareroot(30*one))/(35*one)) );
    const Scalar w3p = as<Scalar>( w2p*(one/(6*one)-ST::squareroot(30*one)/(24*one)) );
    const Scalar w4p = as<Scalar>( w2p*(one/(21*one)-5*one*ST::squareroot(30*one)/(168*one)) );
    const Scalar w5p = as<Scalar>( w2p-2*w3p );
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
    int order = 8;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 4 Stage 8th Order Kuntzmann & Butcher"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage4thOrderHammerHollingsworth_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage4thOrderHammerHollingsworth_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Hammer & Hollingsworth method\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.3, pg 207\n"
                << "c = [ 1/2-sqrt(3)/6  1/2+sqrt(3)/6 ]'\n"
                << "A = [ 1/4            1/4-sqrt(3)/6 ]\n"
                << "    [ 1/4+sqrt(3)/6  1/4           ]\n"
                << "b = [ 1/2            1/2           ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar onequarter = as<Scalar>( one/(4*one) );
    const Scalar onehalf = as<Scalar>( one/(2*one) );
    A(0,0) = onequarter;
    A(0,1) = as<Scalar>( onequarter-ST::squareroot(3*one)/(6*one) );
    A(1,0) = as<Scalar>( onequarter+ST::squareroot(3*one)/(6*one) );
    A(1,1) = onequarter;
    b(0) = onehalf;
    b(1) = onehalf;
    c(0) = as<Scalar>( onehalf - ST::squareroot(3*one)/(6*one) );
    c(1) = as<Scalar>( onehalf + ST::squareroot(3*one)/(6*one) );
    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 4th Order Hammer & Hollingsworth"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class IRK1StageTheta_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  IRK1StageTheta_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Non-standard finite-difference methods\n"
                << "in dynamical systems, P. Kama,\n"
                << "Dissertation, University of Pretoria, pg. 49.\n"
                << "Comment:  Generalized Implicit Midpoint Method\n"
                << "c = [ theta ]'\n"
                << "A = [ theta ]\n"
                << "b = [  1  ]'\n"
                << std::endl;

    typedef Teuchos::ScalarTraits<Scalar> ST;
    theta_default_ = ST::one()/(2*ST::one());
    theta_ = theta_default_;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description() and
      pl->get<std::string>("Stepper Type") != "Implicit Midpoint"
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));
    theta_ = pl->get<double>("theta",theta_default_);

    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    A(0,0) = theta_;
    b(0) = ST::one();
    c(0) = theta_;

    int order = 1;
    if (theta_ == theta_default_) order = 2;

    this->initialize(A, b, c, order, 1, 2, this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  virtual std::string description() const { return "IRK 1 Stage Theta Method"; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");
    pl->set<double>("theta",theta_default_,
      "Valid values are 0 <= theta <= 1, where theta = 0 "
      "implies Forward Euler, theta = 1/2 implies implicit midpoint "
      "method (default), and theta = 1 implies Backward Euler. "
      "For theta != 1/2, this method is first-order accurate, "
      "and with theta = 1/2, it is second-order accurate.  "
      "This method is A-stable, but becomes L-stable with theta=1.");

    return pl;
  }

  private:
    Scalar theta_default_;
    Scalar theta_;
};


// ----------------------------------------------------------------------------
template<class Scalar>
class EDIRK2StageTheta_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  EDIRK2StageTheta_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "Computer Methods for ODEs and DAEs\n"
                << "U. M. Ascher and L. R. Petzold\n"
                << "p. 113\n"
                << "c = [  0       1     ]'\n"
                << "A = [  0       0     ]\n"
                << "    [ 1-theta  theta ]\n"
                << "b = [ 1-theta  theta ]'\n"
                << std::endl;

    typedef Teuchos::ScalarTraits<Scalar> ST;
    theta_default_ = ST::one()/(2*ST::one());
    theta_ = theta_default_;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));
    theta_ = pl->get<double>("theta",theta_default_);

    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    TEUCHOS_TEST_FOR_EXCEPTION(
      theta_ == zero, std::logic_error,
      "'theta' can not be zero, as it makes this IRK stepper explicit.");

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    A(0,0) = zero;
    A(0,1) = zero;
    A(1,0) = Teuchos::as<Scalar>( one - theta_ );
    A(1,1) = theta_;
    b(0) = Teuchos::as<Scalar>( one - theta_ );
    b(1) = theta_;
    c(0) = theta_;
    c(1) = one;

    int order = 1;
    if (theta_ == theta_default_) order = 2;

    this->initialize(A, b, c, order, 1, 2, this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");
    pl->set<double>("theta",theta_default_,
      "Valid values are 0 < theta <= 1, where theta = 0 "
      "implies Forward Euler, theta = 1/2 implies trapezoidal "
      "method (default), and theta = 1 implies Backward Euler. "
      "For theta != 1/2, this method is first-order accurate, "
      "and with theta = 1/2, it is second-order accurate.  "
      "This method is A-stable, but becomes L-stable with theta=1.");

    return pl;
  }

  virtual std::string description() const {return "EDIRK 2 Stage Theta Method";}

  private:
    Scalar theta_default_;
    Scalar theta_;
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit1Stage2ndOrderGauss_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit1Stage2ndOrderGauss_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.2, pg 72\n"
                << "Also:  Implicit midpoint rule\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.1, pg 205\n"
                << "c = [ 1/2 ]'\n"
                << "A = [ 1/2 ]\n"
                << "b = [  1  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar onehalf = ST::one()/(2*ST::one());
    const Scalar one = ST::one();
    A(0,0) = onehalf;
    b(0) = one;
    c(0) = onehalf;
    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 1 Stage 2nd order Gauss"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage4thOrderGauss_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage4thOrderGauss_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.2, pg 72\n"
                << "c = [ 1/2-sqrt(3)/6  1/2+sqrt(3)/6 ]'\n"
                << "A = [ 1/4            1/4-sqrt(3)/6 ]\n"
                << "    [ 1/4+sqrt(3)/6  1/4           ]\n"
                << "b = [ 1/2            1/2 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar onehalf = as<Scalar>(one/(2*one));
    const Scalar three = as<Scalar>(3*one);
    const Scalar six = as<Scalar>(6*one);
    const Scalar onefourth = as<Scalar>(one/(4*one));
    const Scalar alpha = ST::squareroot(three)/six;

    A(0,0) = onefourth;
    A(0,1) = onefourth-alpha;
    A(1,0) = onefourth+alpha;
    A(1,1) = onefourth;
    b(0) = onehalf;
    b(1) = onehalf;
    c(0) = onehalf-alpha;
    c(1) = onehalf+alpha;
    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 4th order Gauss"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage6thOrderGauss_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage6thOrderGauss_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "Table 5.2, pg 72\n"
      << "c = [ 1/2-sqrt(15)/10   1/2              1/2+sqrt(15)/10  ]'\n"
      << "A = [ 5/36              2/9-sqrt(15)/15  5/36-sqrt(15)/30 ]\n"
      << "    [ 5/36+sqrt(15)/24  2/9              5/36-sqrt(15)/24 ]\n"
      << "    [ 5/36+sqrt(15)/30  2/9+sqrt(15)/15  5/36             ]\n"
      << "b = [ 5/18              4/9              5/18             ]'"
      << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar ten = as<Scalar>(10*one);
    const Scalar fifteen = as<Scalar>(15*one);
    const Scalar twentyfour = as<Scalar>(24*one);
    const Scalar thirty = as<Scalar>(30*one);
    const Scalar sqrt15over10 = as<Scalar>(ST::squareroot(fifteen)/ten);
    const Scalar sqrt15over15 = as<Scalar>(ST::squareroot(fifteen)/fifteen);
    const Scalar sqrt15over24 = as<Scalar>(ST::squareroot(fifteen)/twentyfour);
    const Scalar sqrt15over30 = as<Scalar>(ST::squareroot(fifteen)/thirty);

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
    int order = 6;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 6th order Gauss"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit1Stage1stOrderRadauA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit1Stage1stOrderRadauA_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.3, pg 73\n"
                << "c = [ 0 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    A(0,0) = one;
    b(0) = one;
    c(0) = zero;
    int order = 1;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 1 Stage 1st order Radau left"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage3rdOrderRadauA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage3rdOrderRadauA_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.3, pg 73\n"
                << "c = [  0    2/3 ]'\n"
                << "A = [ 1/4  -1/4 ]\n"
                << "    [ 1/4  5/12 ]\n"
                << "b = [ 1/4  3/4  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    A(0,0) = as<Scalar>(one/(4*one));
    A(0,1) = as<Scalar>(-one/(4*one));
    A(1,0) = as<Scalar>(one/(4*one));
    A(1,1) = as<Scalar>(5*one/(12*one));
    b(0) = as<Scalar>(one/(4*one));
    b(1) = as<Scalar>(3*one/(4*one));
    c(0) = zero;
    c(1) = as<Scalar>(2*one/(3*one));
    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 3rd order Radau left"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage5thOrderRadauA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage5thOrderRadauA_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.4, pg 73\n"
                << "c = [  0   (6-sqrt(6))/10       (6+sqrt(6))/10      ]'\n"
                << "A = [ 1/9  (-1-sqrt(6))/18      (-1+sqrt(6))/18     ]\n"
                << "    [ 1/9  (88+7*sqrt(6))/360   (88-43*sqrt(6))/360 ]\n"
                << "    [ 1/9  (88+43*sqrt(6))/360  (88-7*sqrt(6))/360  ]\n"
                << "b = [ 1/9  (16+sqrt(6))/36      (16-sqrt(6))/36     ]'"
                << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 5;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 5th order Radau left"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit1Stage1stOrderRadauB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit1Stage1stOrderRadauB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.5, pg 74\n"
                << "c = [ 1 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    A(0,0) = one;
    b(0) = one;
    c(0) = one;
    int order = 1;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 1 Stage 1st order Radau right"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage3rdOrderRadauB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage3rdOrderRadauB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.5, pg 74\n"
                << "c = [ 1/3     1   ]'\n"
                << "A = [ 5/12  -1/12 ]\n"
                << "    [ 3/4    1/4  ]\n"
                << "b = [ 3/4    1/4  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    A(0,0) = as<Scalar>( 5*one/(12*one) );
    A(0,1) = as<Scalar>( -one/(12*one) );
    A(1,0) = as<Scalar>( 3*one/(4*one) );
    A(1,1) = as<Scalar>( one/(4*one) );
    b(0) = as<Scalar>( 3*one/(4*one) );
    b(1) = as<Scalar>( one/(4*one) );
    c(0) = as<Scalar>( one/(3*one) );
    c(1) = one;
    int order = 3;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 3rd order Radau right"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage5thOrderRadauB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage5thOrderRadauB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "Table 5.6, pg 74\n"
      << "c = [ (4-sqrt(6))/10          (4+sqrt(6))/10          1    ]'\n"
      << "A = [ A1 A2 A3 ]\n"
      << "      A1 = [ (88-7*sqrt(6))/360     ]\n"
      << "           [ (296+169*sqrt(6))/1800 ]\n"
      << "           [ (16-sqrt(6))/36        ]\n"
      << "      A2 = [ (296-169*sqrt(6))/1800 ]\n"
      << "           [ (88+7*sqrt(6))/360     ]\n"
      << "           [ (16+sqrt(6))/36        ]\n"
      << "      A3 = [ (-2+3*sqrt(6))/225 ]\n"
      << "           [ (-2-3*sqrt(6))/225 ]\n"
      << "           [ 1/9                ]\n"
      << "b = [ (16-sqrt(6))/36         (16+sqrt(6))/36         1/9 ]'"
      << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
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
    int order = 5;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 5th order Radau right"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage2ndOrderLobattoA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage2ndOrderLobattoA_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.7, pg 75\n"
                << "c = [  0    1   ]'\n"
                << "A = [  0    0   ]\n"
                << "    [ 1/2  1/2  ]\n"
                << "b = [ 1/2  1/2  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    A(0,0) = zero;
    A(0,1) = zero;
    A(1,0) = as<Scalar>( one/(2*one) );
    A(1,1) = as<Scalar>( one/(2*one) );
    b(0) = as<Scalar>( one/(2*one) );
    b(1) = as<Scalar>( one/(2*one) );
    c(0) = zero;
    c(1) = one;
    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 2nd order Lobatto A"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage4thOrderLobattoA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage4thOrderLobattoA_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.7, pg 75\n"
                << "c = [  0    1/2    1  ]'\n"
                << "A = [  0     0     0   ]\n"
                << "    [ 5/24  1/3  -1/24  ]\n"
                << "    [ 1/6   2/3   1/6   ]\n"
                << "b = [ 1/6   2/3   1/6   ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 4th order Lobatto A"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit4Stage6thOrderLobattoA_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit4Stage6thOrderLobattoA_RKBT()
  {
    using Teuchos::as;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "Table 5.8, pg 75\n"
      << "c = [ 0  (5-sqrt(5))/10  (5+sqrt(5))/10  1 ]'\n"
      << "A = [ A1  A2  A3  A4 ]\n"
      << "      A1 = [ 0               ]\n"
      << "           [ (11+sqrt(5)/120 ]\n"
      << "           [ (11-sqrt(5)/120 ]\n"
      << "           [ 1/12            ]\n"
      << "      A2 = [ 0                    ]\n"
      << "           [ (25-sqrt(5))/120     ]\n"
      << "           [ (25+13*sqrt(5))/120  ]\n"
      << "           [ 5/12                 ]\n"
      << "      A3 = [ 0                   ]\n"
      << "           [ (25-13*sqrt(5))/120 ]\n"
      << "           [ (25+sqrt(5))/120    ]\n"
      << "           [ 5/12                ]\n"
      << "      A4 = [ 0                ]\n"
      << "           [ (-1+sqrt(5))/120 ]\n"
      << "           [ (-1-sqrt(5))/120 ]\n"
      << "           [ 1/12             ]\n"
      << "b = [ 1/12  5/12  5/12   1/12 ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 6;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 4 Stage 6th order Lobatto A"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage2ndOrderLobattoB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage2ndOrderLobattoB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.9, pg 76\n"
                << "c = [  0    1   ]'\n"
                << "A = [ 1/2   0   ]\n"
                << "    [ 1/2   0   ]\n"
                << "b = [ 1/2  1/2  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    A(0,0) = as<Scalar>( one/(2*one) );
    A(0,1) = zero;
    A(1,0) = as<Scalar>( one/(2*one) );
    A(1,1) = zero;
    b(0) = as<Scalar>( one/(2*one) );
    b(1) = as<Scalar>( one/(2*one) );
    c(0) = zero;
    c(1) = one;
    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 2nd order Lobatto B"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage4thOrderLobattoB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage4thOrderLobattoB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.9, pg 76\n"
                << "c = [  0    1/2    1   ]'\n"
                << "A = [ 1/6  -1/6    0   ]\n"
                << "    [ 1/6   1/3    0   ]\n"
                << "    [ 1/6   5/6    0   ]\n"
                << "b = [ 1/6   2/3   1/6  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 4th order Lobatto B"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit4Stage6thOrderLobattoB_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit4Stage6thOrderLobattoB_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "Table 5.10, pg 76\n"
      << "c = [ 0     (5-sqrt(5))/10       (5+sqrt(5))/10       1     ]'\n"
      << "A = [ 1/12  (-1-sqrt(5))/24      (-1+sqrt(5))/24      0     ]\n"
      << "    [ 1/12  (25+sqrt(5))/120     (25-13*sqrt(5))/120  0     ]\n"
      << "    [ 1/12  (25+13*sqrt(5))/120  (25-sqrt(5))/120     0     ]\n"
      << "    [ 1/12  (11-sqrt(5))/24      (11+sqrt(5))/24      0     ]\n"
      << "b = [ 1/12  5/12                 5/12                 1/12  ]'"
      << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    A(0,0) = as<Scalar>( one/(12*one) );
    A(0,1) = as<Scalar>( (-one-ST::squareroot(5*one))/(24*one) );
    A(0,2) = as<Scalar>( (-one+ST::squareroot(5*one))/(24*one) );
    A(0,3) = zero;
    A(1,0) = as<Scalar>( one/(12*one) );
    A(1,1) = as<Scalar>( (25*one+ST::squareroot(5*one))/(120*one) );
    A(1,2) = as<Scalar>( (25*one-13*one*ST::squareroot(5*one))/(120*one) );
    A(1,3) = zero;
    A(2,0) = as<Scalar>( one/(12*one) );
    A(2,1) = as<Scalar>( (25*one+13*one*ST::squareroot(5*one))/(120*one) );
    A(2,2) = as<Scalar>( (25*one-ST::squareroot(5*one))/(120*one) );
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
    int order = 6;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 4 Stage 6th order Lobatto B"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit2Stage2ndOrderLobattoC_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit2Stage2ndOrderLobattoC_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.11, pg 76\n"
                << "c = [  0    1   ]'\n"
                << "A = [ 1/2 -1/2  ]\n"
                << "    [ 1/2  1/2  ]\n"
                << "b = [ 1/2  1/2  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    A(0,0) = as<Scalar>( one/(2*one) );
    A(0,1) = as<Scalar>( -one/(2*one) );
    A(1,0) = as<Scalar>( one/(2*one) );
    A(1,1) = as<Scalar>( one/(2*one) );
    b(0) = as<Scalar>( one/(2*one) );
    b(1) = as<Scalar>( one/(2*one) );
    c(0) = zero;
    c(1) = one;
    int order = 2;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 2 Stage 2nd order Lobatto C"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit3Stage4thOrderLobattoC_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit3Stage4thOrderLobattoC_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.11, pg 76\n"
                << "c = [  0    1/2    1   ]'\n"
                << "A = [ 1/6  -1/3   1/6  ]\n"
                << "    [ 1/6   5/12 -1/12 ]\n"
                << "    [ 1/6   2/3   1/6  ]\n"
                << "b = [ 1/6   2/3   1/6  ]'" << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 4;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 3 Stage 4th order Lobatto C"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class Implicit4Stage6thOrderLobattoC_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  Implicit4Stage6thOrderLobattoC_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "Table 5.12, pg 76\n"
      << "c = [ 0     (5-sqrt(5))/10       (5+sqrt(5))/10       1          ]'\n"
      << "A = [ 1/12  -sqrt(5)/12          sqrt(5)/12          -1/12       ]\n"
      << "    [ 1/12  1/4                  (10-7*sqrt(5))/60   sqrt(5)/60  ]\n"
      << "    [ 1/12  (10+7*sqrt(5))/60    1/4                 -sqrt(5)/60 ]\n"
      << "    [ 1/12  5/12                 5/12                 1/12       ]\n"
      << "b = [ 1/12  5/12                 5/12                 1/12       ]'"
      << std::endl;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
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
    int order = 6;

    this->initialize(A,b,c,order,Description.str());
  }
  virtual std::string description() const
    { return "RK Implicit 4 Stage 6th order Lobatto C"; }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class SDIRK5Stage5thOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK5Stage5thOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "A-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "pg101 \n"
      << "c = [ (6-sqrt(6))/10   ]\n"
      << "    [ (6+9*sqrt(6))/35 ]\n"
      << "    [ 1                ]\n"
      << "    [ (4-sqrt(6))/10   ]\n"
      << "    [ (4+sqrt(6))/10   ]\n"
      << "A = [ A1 A2 A3 A4 A5 ]\n"
      << "      A1 = [ (6-sqrt(6))/10               ]\n"
      << "           [ (-6+5*sqrt(6))/14            ]\n"
      << "           [ (888+607*sqrt(6))/2850       ]\n"
      << "           [ (3153-3082*sqrt(6))/14250    ]\n"
      << "           [ (-32583+14638*sqrt(6))/71250 ]\n"
      << "      A2 = [ 0                           ]\n"
      << "           [ (6-sqrt(6))/10              ]\n"
      << "           [ (126-161*sqrt(6))/1425      ]\n"
      << "           [ (3213+1148*sqrt(6))/28500   ]\n"
      << "           [ (-17199+364*sqrt(6))/142500 ]\n"
      << "      A3 = [ 0                       ]\n"
      << "           [ 0                       ]\n"
      << "           [ (6-sqrt(6))/10          ]\n"
      << "           [ (-267+88*sqrt(6))/500   ]\n"
      << "           [ (1329-544*sqrt(6))/2500 ]\n"
      << "      A4 = [ 0                     ]\n"
      << "           [ 0                     ]\n"
      << "           [ 0                     ]\n"
      << "           [ (6-sqrt(6))/10        ]\n"
      << "           [ (-96+131*sqrt(6))/625 ]\n"
      << "      A5 = [ 0              ]\n"
      << "           [ 0              ]\n"
      << "           [ 0              ]\n"
      << "           [ 0              ]\n"
      << "           [ (6-sqrt(6))/10 ]\n"
      << "b = [               0 ]\n"
      << "    [               0 ]\n"
      << "    [             1/9 ]\n"
      << "    [ (16-sqrt(6))/36 ]\n"
      << "    [ (16+sqrt(6))/36 ]" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 5 Stage 5th order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    const Scalar sqrt6 = ST::squareroot(as<Scalar>(6*one));
    const Scalar gamma = as<Scalar>( (6*one - sqrt6) / (10*one) ); // diagonal
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

    int order = 5;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class SDIRK5Stage4thOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK5Stage4thOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
      << "L-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
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

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 5 Stage 4th order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    const Scalar onequarter = as<Scalar>( one/(4*one) );
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

    int order = 4;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};


// ----------------------------------------------------------------------------
template<class Scalar>
class SDIRK3Stage4thOrder_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK3Stage4thOrder_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "pg100 \n"
                << "gamma = (1/sqrt(3))*cos(pi/18)+1/2\n"
                << "delta = 1/(6*(2*gamma-1)^2)\n"
                << "c = [ gamma      1/2        1-gamma ]'\n"
                << "A = [ gamma                         ]\n"
                << "    [ 1/2-gamma  gamma              ]\n"
                << "    [ 2*gamma    1-4*gamma  gamma   ]\n"
                << "b = [ delta      1-2*delta  delta   ]'" << std::endl;

    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 3 Stage 4th order"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    const Scalar pi = as<Scalar>(4*one)*std::atan(one);
    const Scalar gamma = as<Scalar>( one/ST::squareroot(3*one)*std::cos(pi/(18*one))+one/(2*one) );
    const Scalar delta = as<Scalar>( one/(6*one*std::pow(2*gamma-one,2*one)) );
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

    int order = 4;

    this->initialize(A,b,c,order,this->getDescription());
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};

// ----------------------------------------------------------------------------
/** \brief SDIRK 2(1) pair
 *
 *  The tableau (order=2(1)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\ \hline
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0 & 0   & \\
 *                         1 & -1  & 1 \\ \hline
 *                           & 1/2 & 1/2 \\
 *                           & 1   & 0 \end{array}
 *  \f]
 *
 */
template<class Scalar>
class SDIRK21_RKBT :
  virtual public RKButcherTableau<Scalar>
{
  public:
  SDIRK21_RKBT()
  {
    std::ostringstream Description;
    Description << this->description() << "\n"
                << "c =     [  1  0   ]'\n"
                << "A =     [  1      ]\n"
                << "        [ -1  1   ]\n"
                << "b     = [ 1/2 1/2 ]\n"
                << "bstar = [  1  0   ]\n" << std::endl;
    this->setDescription(Description.str());
    this->setParameterList(Teuchos::null);
  }

  virtual std::string description() const { return "SDIRK 2(1) Pair"; }

  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& pList)
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    if (pList == Teuchos::null) *pl = *(this->getValidParameters());
    else pl = pList;
    // Can not validate because optional parameters (e.g., Solver Name).
    //pl->validateParametersAndSetDefaults(*this->getValidParameters());
    TEUCHOS_TEST_FOR_EXCEPTION(
      pl->get<std::string>("Stepper Type") != this->description()
      ,std::runtime_error,
      "  Stepper Type != \""+this->description()+"\"\n"
      "  Stepper Type = " + pl->get<std::string>("Stepper Type"));

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) = one;
    A(0,1) = zero;

    //A(1,0) =
    A(1,0) = -one;
    A(1,1) =  one;

    // Fill b:
    b(0) = as<Scalar>(one/(2*one));
    b(1) = as<Scalar>(one/(2*one));

    // Fill c:
    c(0) = one;
    c(1) = zero;

    // Fill bstar
    bstar(0) = one;
    bstar(1) = zero;
    int order = 2;

    this->initialize(A,b,c,order,this->getDescription(),true,bstar);
    this->setMyParamList(pl);
    this->rkbtPL_ = pl;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    pl->setName("Default Stepper - " + this->description());
    pl->set<std::string>("Description", this->getDescription());
    pl->set<std::string>("Stepper Type", this->description());
    pl->set<bool>("Use Embedded", false);
    pl->set<std::string>("Solver Name", "",
      "Name of ParameterList containing the solver specifications.");

    return pl;
  }
};


} // namespace Tempus


#endif // Tempus_RKButcherTableau_hpp
