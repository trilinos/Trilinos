// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_StepperRKButcherTableau_hpp
#define Tempus_StepperRKButcherTableau_hpp

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include "Tempus_StepperExplicitRK.hpp"
#include "Tempus_StepperDIRK.hpp"
#include "Tempus_RKButcherTableau.hpp"


namespace Tempus {


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
 *  \begin{array}{c|c} 0 & 0 \\ \hline
 *                       & 1 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERK_ForwardEuler :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_ForwardEuler()
  {
    this->setStepperType("RK Forward Euler");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_ForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Forward Euler");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_ForwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Forward Euler");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c = [ 0 ]'\n"
                << "A = [ 0 ]\n"
                << "b = [ 1 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::SerialDenseMatrix<int,Scalar> A(1,1);
    Teuchos::SerialDenseVector<int,Scalar> b(1);
    Teuchos::SerialDenseVector<int,Scalar> c(1);
    A(0,0) = ST::zero();
    b(0) = ST::one();
    c(0) = ST::zero();
    int order = 1;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
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
class StepperERK_4Stage4thOrder :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_4Stage4thOrder()
  {
    this->setStepperType("RK Explicit 4 Stage");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_4Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 4 Stage");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_4Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 4 Stage");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
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
                << "b = [ 1/6 1/3 1/3 1/6 ]'";
    return Description.str();
  }

  void setupTableau()
  {
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
    A(0,0) =    zero; A(0,1) =    zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) = onehalf; A(1,1) =    zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) =    zero; A(2,1) = onehalf; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =    zero; A(3,1) =    zero; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = onethird; b(2) = onethird; b(3) = onesixth;

    // fill c:
    c(0) = zero; c(1) = onehalf; c(2) = onehalf; c(3) = one;

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RK Bogacki-Shampine Butcher Tableau
 *
 *  The tableau (order=3(2)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A   \\ \hline
 *      & b^T \\
 *      & b^{*T}
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0  & 0    &     &     & \\
 *                        1/2 & 1/2  & 0   &     & \\
 *                        3/4 & 0    & 3/4 & 0   & \\
 *                         1  & 2/9  & 1/3 & 4/9 & 0 \\ \hline
 *                            & 2/9  & 1/3 & 4/9 & 0 \\
 *                            & 7/24 & 1/4 & 1/3 & 1/8 \end{array}
 *  \f]
 *  Reference:  P. Bogacki and L.F. Shampine.
 *              A 3(2) pair of Runge–Kutta formulas.
 *              Applied Mathematics Letters, 2(4):321 – 325, 1989.
 */
template<class Scalar>
class StepperERK_BogackiShampine32 :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_BogackiShampine32()
  {
    this->setStepperType("Bogacki-Shampine 3(2) Pair");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_BogackiShampine32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("Bogacki-Shampine 3(2) Pair");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_BogackiShampine32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("Bogacki-Shampine 3(2) Pair");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "P. Bogacki and L.F. Shampine.\n"
                << "A 3(2) pair of Runge–Kutta formulas.\n"
                << "Applied Mathematics Letters, 2(4):321 – 325, 1989.\n"
                << "c =     [ 0     1/2  3/4   1  ]'\n"
                << "A =     [ 0                   ]\n"
                << "        [ 1/2    0            ]\n"
                << "        [  0    3/4   0       ]\n"
                << "        [ 2/9   1/3  4/9   0  ]\n"
                << "b     = [ 2/9   1/3  4/9   0  ]'\n"
                << "bstar = [ 7/24  1/4  1/3  1/8 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);
    const Scalar onethird = one/(3*one);
    const Scalar threefourths = (3*one)/(4*one);
    const Scalar twoninths = (2*one)/(9*one);
    const Scalar fourninths = (4*one)/(9*one);

    // Fill A:
    A(0,0) =     zero; A(0,1) =        zero; A(0,2) =      zero; A(0,3) = zero;
    A(1,0) =  onehalf; A(1,1) =        zero; A(1,2) =      zero; A(1,3) = zero;
    A(2,0) =     zero; A(2,1) =threefourths; A(2,2) =      zero; A(2,3) = zero;
    A(3,0) =twoninths; A(3,1) =    onethird; A(3,2) =fourninths; A(3,3) = zero;

    // Fill b:
    b(0) = A(3,0); b(1) = A(3,1); b(2) = A(3,2); b(3) = A(3,3);

    // Fill c:
    c(0) = zero; c(1) = onehalf; c(2) = threefourths; c(3) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(7*one/(24*one));
    bstar(1) = as<Scalar>(1*one/(4*one));
    bstar(2) = as<Scalar>(1*one/(3*one));
    bstar(3) = as<Scalar>(1*one/(8*one));
    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
};


// ----------------------------------------------------------------------------
/** \brief Explicit RK Merson Butcher Tableau
 *
 *  The tableau (order=4(5)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A   \\ \hline
 *      & b^T \\
 *      & b^{*T}
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
class StepperERK_Merson45 :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_Merson45()
  {
    this->setStepperType("Merson 4(5) Pair");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_Merson45(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("Merson 4(5) Pair");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_Merson45(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("Merson 4(5) Pair");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
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
                << "b     = [ 1/6    0    0   2/3  1/6 ]'\n"
                << "bstar = [ 1/10   0  3/10  2/5  1/5 ]'";
    return Description.str();
  }


protected:

  void setupTableau()
  {
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

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
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
class StepperERK_3_8Rule :
  virtual public StepperExplicitRK<Scalar>
{
public:

  StepperERK_3_8Rule()
  {
    this->setStepperType("RK Explicit 3/8 Rule");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_3_8Rule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 3/8 Rule");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_3_8Rule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 3/8 Rule");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.2, pg 138\n"
                << "c = [  0  1/3 2/3  1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/3  0          ]\n"
                << "    [-1/3  1   0      ]\n"
                << "    [  1  -1   1   0  ]\n"
                << "b = [ 1/8 3/8 3/8 1/8 ]'";
    return Description.str();
  }


protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onethird     = as<Scalar>(one/(3*one));
    const Scalar twothirds    = as<Scalar>(2*one/(3*one));
    const Scalar oneeighth    = as<Scalar>(one/(8*one));
    const Scalar threeeighths = as<Scalar>(3*one/(8*one));

    // Fill A:
    A(0,0) =      zero; A(0,1) = zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) =  onethird; A(1,1) = zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) = -onethird; A(2,1) =  one; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =       one; A(3,1) = -one; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) =oneeighth; b(1) =threeeighths; b(2) =threeeighths; b(3) =oneeighth;

    // Fill c:
    c(0) = zero; c(1) = onethird; c(2) = twothirds; c(3) = one;

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
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
class StepperERK_4Stage3rdOrderRunge :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_4Stage3rdOrderRunge()
  {
    this->setStepperType("RK Explicit 4 Stage 3rd order by Runge");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_4Stage3rdOrderRunge(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 4 Stage 3rd order by Runge");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_4Stage3rdOrderRunge(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 4 Stage 3rd order by Runge");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2  1   1  ]'\n"
                << "A = [  0              ]\n"
                << "    [ 1/2  0          ]\n"
                << "    [  0   1   0      ]\n"
                << "    [  0   0   1   0  ]\n"
                << "b = [ 1/6 2/3  0  1/6 ]'";
    return Description.str();
  }
protected:

  void setupTableau()
  {
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
    A(0,0) =    zero; A(0,1) = zero; A(0,2) = zero; A(0,3) = zero;
    A(1,0) = onehalf; A(1,1) = zero; A(1,2) = zero; A(1,3) = zero;
    A(2,0) =    zero; A(2,1) =  one; A(2,2) = zero; A(2,3) = zero;
    A(3,0) =    zero; A(3,1) = zero; A(3,2) =  one; A(3,3) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = twothirds; b(2) = zero; b(3) = onesixth;

    // Fill c:
    c(0) = zero; c(1) = onehalf; c(2) = one; c(3) = one;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
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
class StepperERK_5Stage3rdOrderKandG :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_5Stage3rdOrderKandG()
  {
    this->setStepperType("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_5Stage3rdOrderKandG(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_5Stage3rdOrderKandG(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 5 Stage 3rd order by Kinnmark and Gray");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Kinnmark & Gray 5 stage, 3rd order scheme \n"
                << "Modified by P. Ullrich.  From the prim_advance_mod.F90 \n"
                << "routine in the HOMME atmosphere model code.\n"
                << "c = [  0  1/5  1/5  1/3  2/3  ]'\n"
                << "A = [  0                      ]\n"
                << "    [ 1/5  0                  ]\n"
                << "    [  0  1/5   0             ]\n"
                << "    [  0   0   1/3   0        ]\n"
                << "    [  0   0    0   2/3   0   ]\n"
                << "b = [ 1/4  0    0    0   3/4  ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
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
    A(0,0) =     zero; A(0,1) =     zero; A(0,2) =     zero; A(0,3) =      zero; A(0,4) = zero;
    A(1,0) = onefifth; A(1,1) =     zero; A(1,2) =     zero; A(1,3) =      zero; A(1,4) = zero;
    A(2,0) =     zero; A(2,1) = onefifth; A(2,2) =     zero; A(2,3) =      zero; A(2,4) = zero;
    A(3,0) =     zero; A(3,1) =     zero; A(3,2) = onethird; A(3,3) =      zero; A(3,4) = zero;
    A(4,0) =     zero; A(4,1) =     zero; A(4,2) =     zero; A(4,3) = twothirds; A(4,4) = zero;

    // Fill b:
    b(0) =onefourth; b(1) =zero; b(2) =zero; b(3) =zero; b(4) =threefourths;

    // Fill c:
    c(0) =zero; c(1) =onefifth; c(2) =onefifth; c(3) =onethird; c(4) =twothirds;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
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
class StepperERK_3Stage3rdOrder :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_3Stage3rdOrder()
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_3Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_3Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c = [  0  1/2  1  ]'\n"
                << "A = [  0          ]\n"
                << "    [ 1/2  0      ]\n"
                << "    [ -1   2   0  ]\n"
                << "b = [ 1/6 4/6 1/6 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
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
    A(0,0) =    zero; A(0,1) = zero; A(0,2) = zero;
    A(1,0) = onehalf; A(1,1) = zero; A(1,2) = zero;
    A(2,0) =    -one; A(2,1) =  two; A(2,2) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = foursixth; b(2) = onesixth;

    // fill c:
    c(0) = zero; c(1) = onehalf; c(2) = one;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
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
class StepperERK_3Stage3rdOrderTVD :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_3Stage3rdOrderTVD()
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order TVD");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_3Stage3rdOrderTVD(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order TVD");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_3Stage3rdOrderTVD(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order TVD");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                  << "This Stepper is known as 'RK Explicit 3 Stage 3rd order TVD' or 'SSPERK33'.\n"
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
                  << "u^(n+1) = u^n/3 + 2 u2/2 + 2 dt L(u2)/3";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
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
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    // Fill A:
    A(0,0) =      zero; A(0,1) =      zero; A(0,2) = zero;
    A(1,0) =       one; A(1,1) =      zero; A(1,2) = zero;
    A(2,0) = onefourth; A(2,1) = onefourth; A(2,2) = zero;

    // Fill b:
    b(0) = onesixth; b(1) = onesixth; b(2) = foursixth;

    // fill c:
    c(0) = zero; c(1) = one; c(2) = onehalf;

    // Fill bstar:
    bstar(0) = as<Scalar>(0.291485418878409);
    bstar(1) = as<Scalar>(0.291485418878409);
    bstar(2) = as<Scalar>(0.417029162243181);

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
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
class StepperERK_3Stage3rdOrderHeun :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_3Stage3rdOrderHeun()
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order by Heun");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_3Stage3rdOrderHeun(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order by Heun");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_3Stage3rdOrderHeun(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit 3 Stage 3rd order by Heun");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/3 2/3 ]'\n"
                << "A = [  0          ] \n"
                << "    [ 1/3  0      ]\n"
                << "    [  0  2/3  0  ]\n"
                << "b = [ 1/4  0  3/4 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
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
    A(0,0) =     zero; A(0,1) =      zero; A(0,2) = zero;
    A(1,0) = onethird; A(1,1) =      zero; A(1,2) = zero;
    A(2,0) =     zero; A(2,1) = twothirds; A(2,2) = zero;

    // Fill b:
    b(0) = onefourth; b(1) = zero; b(2) = threefourths;

    // fill c:
    c(0) = zero; c(1) = onethird; c(2) = twothirds;

    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit Midpoint
 *
 *  The tableau (order=2) is
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
class StepperERK_Midpoint :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_Midpoint()
  {
    this->setStepperType("RK Explicit Midpoint");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_Midpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit Midpoint");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_Midpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit Midpoint");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S.P. Norsett, G. Wanner\n"
                << "Table 1.1, pg 135\n"
                << "c = [  0  1/2 ]'\n"
                << "A = [  0      ]\n"
                << "    [ 1/2  0  ]\n"
                << "b = [  0   1  ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =    zero; A(0,1) = zero;
    A(1,0) = onehalf; A(1,1) = zero;

    // Fill b:
    b(0) = zero; b(1) = one;

    // fill c:
    c(0) = zero; c(1) = onehalf;

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit Ralston
 *
 *  The tableau (order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}   0   &   0   &     \\
 *                       2/3  &  2/3  &  0  \\ \hline
 *                            &  1/4  & 3/4  \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERK_Ralston :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_Ralston()
  {
    this->setStepperType("RK Explicit Ralston");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_Ralston(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit Ralston");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_Ralston(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit Ralston");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "This Stepper is known as 'RK Explicit Ralston' or 'RK2'.\n"
                << "c = [   0   2/3 ]'\n"
                << "A = [   0       ]\n"
                << "    [  2/3   0  ]\n"
                << "b = [  1/4  3/4 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
   typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    const int NumStages = 2;
    const int order = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) = zero; A(0,1) = zero; A(1,1) = zero;
    A(1,0) = (2*one)/(3*one);

    // Fill b:
    b(0) = (1*one)/(4*one);
    b(1) = (3*one)/(4*one);

    // fill c:
    c(0) = zero;
    c(1) = (2*one)/(3*one);


    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RK Explicit Trapezoidal
 *
 *  The tableau (order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                       1  &  1  &  0  \\ \hline
 *                          & 1/2 & 1/2 \\
 *                          & 3/4 & 1/4 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperERK_Trapezoidal :
  virtual public StepperExplicitRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperERK_Trapezoidal()
  {
    this->setStepperType("RK Explicit Trapezoidal");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_Trapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("RK Explicit Trapezoidal");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_Trapezoidal(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Explicit Trapezoidal");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "This Stepper is known as 'RK Explicit Trapezoidal' or 'Heuns Method' or 'SSPERK22'.\n"
                << "c      = [  0    1  ]'\n"
                << "A      = [  0       ]\n"
                << "         [  1    0  ]\n"
                << "b      = [ 1/2  1/2 ]\n"
                << "bstart = [ 3/4  1/4 ]'";
    return Description.str();
  }

protected:

  void setupTableau()
  {
   typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = one/(2*one);

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);

    // Fill A:
    A(0,0) = zero; A(0,1) = zero;
    A(1,0) =  one; A(1,1) = zero;

    // Fill b:
    b(0) = onehalf; b(1) = onehalf;

    // fill c:
    c(0) = zero; c(1) = one;

    // Fill bstar
    bstar(0) = as<Scalar>(3*one/(4*one));
    bstar(1) = as<Scalar>(1*one/(4*one));

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Explicit RK Butcher Tableau
 *
 *  The tableau (stage=5, order=4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\
 *      & \hat{b}^T
 *  \end{array}
 *
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 *
 *
 */
template<class Scalar>
class StepperERK_SSPERK54 :
  virtual public StepperExplicitRK<Scalar>
{
  public:
  StepperERK_SSPERK54()
  {
    this->setStepperType("SSPERK54");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_SSPERK54(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded)
  {
    this->setStepperType("SSPERK54");
    this->setupTableau();
    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_SSPERK54(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SSPERK54");
    this->setupTableau();
    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Strong Stability Preserving Explicit RK (stage=5, order=4)\n"
                << std::endl;
    return Description.str();
  }

protected:

  void setupTableau()
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 5;
    const int order     = 4;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> bstar(NumStages);
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) = A(0,1) =  A(0,2) = A(0,3) = A(0,4) = zero;

    A(1,0) = as<Scalar>(0.391752226571889);
    A(1,1) = A(1,2) = A(1,3) = A(0,4) = zero;

    A(2,0) = as<Scalar>(0.217669096261169);
    A(2,1) = as<Scalar>(0.368410593050372);
    A(2,2) = A(2,3) = A(2,4) = zero;

    A(3,0) = as<Scalar>(0.082692086657811);
    A(3,1) = as<Scalar>(0.139958502191896);
    A(3,2) = as<Scalar>(0.251891774271693);
    A(3,3) = A(2,4) = zero;

    A(4,0) = as<Scalar>(0.067966283637115);
    A(4,1) = as<Scalar>(0.115034698504632);
    A(4,2) = as<Scalar>(0.207034898597385);
    A(4,3) = as<Scalar>(0.544974750228520);
    A(4,4) = zero;

    // Fill b:
    b(0) = as<Scalar>(0.146811876084786);
    b(1) = as<Scalar>(0.248482909444976);
    b(2) = as<Scalar>(0.104258830331980);
    b(3) = as<Scalar>(0.274438900901350);
    b(4) = as<Scalar>(0.226007483236908);

    // fill c:
    c(0) = zero;
    c(1) = A(1,0);
    c(2) = A(2,0) + A(2,1);
    c(3) = A(3,0) + A(3,1) + A(3,2);
    c(4) = A(4,0) + A(4,1) + A(4,2) + A(4,3);

    // Fill bstar:
    bstar(0) = as<Scalar>(0.130649104813131);
    bstar(1) = as<Scalar>(0.317716031201302);
    bstar(2) = as<Scalar>(0.000000869337261);
    bstar(3) = as<Scalar>(0.304581512634772);
    bstar(4) = as<Scalar>(0.247052482013534);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
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
class StepperERK_General :
  virtual public StepperExplicitRK<Scalar>
{
public:
  StepperERK_General()
  {
    this->setStepperType("General ERK");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperERK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    const int order,
    const int orderMin,
    const int orderMax,
    const Teuchos::SerialDenseVector<int,Scalar>& bstar)
  {
    this->setStepperType("General ERK");
    this->setTableau(A,b,c,order,orderMin,orderMax,bstar);

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->tableau_->isImplicit() == true, std::logic_error,
      "Error - General ERK received an implicit Butcher Tableau!\n");

    this->setup(appModel, obs, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded);
  }
#endif
  StepperERK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    const int order,
    const int orderMin,
    const int orderMax,
    const Teuchos::SerialDenseVector<int,Scalar>& bstar,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("General ERK");
    this->setTableau(A,b,c,order,orderMin,orderMax,bstar);

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->tableau_->isImplicit() == true, std::logic_error,
      "Error - General ERK received an implicit Butcher Tableau!\n");

    this->setup(appModel, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, stepperRKAppAction);
  }

  virtual std::string getDescription() const
  {
    std::stringstream Description;
    Description << this->getStepperType() << "\n"
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
      << "b = [ 1/6 1/3 1/3 1/6 ]'";
    return Description.str();
  }

  void setupTableau()
  {
    if (this->tableau_ == Teuchos::null) {
      // Set tableau to the default if null, otherwise keep current tableau.
      auto stepper = Teuchos::rcp(new StepperERK_4Stage4thOrder<Scalar>());
      auto t = stepper->getTableau();
      this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
                                 this->getStepperType(),
                                 t->A(),t->b(),t->c(),
                                 t->order(),t->orderMin(),t->orderMax(),
                                 t->bstar()));
    }
  }

  void setTableau(const Teuchos::SerialDenseMatrix<int,Scalar>& A,
                  const Teuchos::SerialDenseVector<int,Scalar>& b,
                  const Teuchos::SerialDenseVector<int,Scalar>& c,
                  const int order,
                  const int orderMin,
                  const int orderMax,
                  const Teuchos::SerialDenseVector<int,Scalar>&
                    bstar = Teuchos::SerialDenseVector<int,Scalar>())
  {
    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,orderMin,orderMax,bstar));
    this->isInitialized_ = false;
  }

  virtual std::string getDefaultICConsistency() const { return "Consistent"; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicERK(pl);
    pl->set<std::string>("Initial Condition Consistency",
                         this->getDefaultICConsistency());

    // Tableau ParameterList
    Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
    tableauPL->set<std::string>("A",
     "0.0 0.0 0.0 0.0; 0.5 0.0 0.0 0.0; 0.0 0.5 0.0 0.0; 0.0 0.0 1.0 0.0");
    tableauPL->set<std::string>("b",
     "0.166666666666667 0.333333333333333 0.333333333333333 0.166666666666667");
    tableauPL->set<std::string>("c", "0.0 0.5 0.5 1.0");
    tableauPL->set<int>("order", 4);
    tableauPL->set<std::string>("bstar", "");
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
class StepperDIRK_BackwardEuler :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   *  Requires subsequent setModel() and initialize()
   *  calls before calling takeStep().
  */
  StepperDIRK_BackwardEuler()
  {
    this->setStepperType("RK Backward Euler");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperDIRK_BackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("RK Backward Euler");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperDIRK_BackwardEuler(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Backward Euler");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c = [ 1 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) = ST::one();

    // Fill b:
    b(0) = ST::one();

    // Fill c:
    c(0) = ST::one();

    int order = 1;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
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
 *  The default value is \f$\gamma = (2 - \sqrt{2})/2\f$.
 *  This will produce an L-stable 2nd order method with the stage
 *  times within the timestep.  Other values of gamma will still
 *  produce an L-stable scheme, but will only be 1st order accurate.
 *  L-stability is guaranteed because \f$A_{sj} = b_j\f$.
 *
 *  Reference: U. M. Ascher and L. R. Petzold,
 *             Computer Methods for ODEs and DAEs, p. 106.
 */
template<class Scalar>
class StepperSDIRK_2Stage2ndOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_2Stage2ndOrder()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    gammaDefault_ = Teuchos::as<Scalar>((2*one-ST::squareroot(2*one))/(2*one));

    this->setStepperType("SDIRK 2 Stage 2nd order");
    this->setGamma(gammaDefault_);
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_2Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    Scalar gamma = Scalar(0.2928932188134524))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    gammaDefault_ = Teuchos::as<Scalar>((2*one-ST::squareroot(2*one))/(2*one));

    this->setStepperType("SDIRK 2 Stage 2nd order");
    this->setGamma(gamma);
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_2Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    Scalar gamma = Scalar(0.2928932188134524))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    gammaDefault_ = Teuchos::as<Scalar>((2*one-ST::squareroot(2*one))/(2*one));

    this->setStepperType("SDIRK 2 Stage 2nd order");
    this->setGamma(gamma);
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  void setGamma(Scalar gamma)
  {
    gamma_ = gamma;
    this->isInitialized_ = false;
    this->setupTableau();
  }

  Scalar getGamma() { return gamma_; }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Computer Methods for ODEs and DAEs\n"
                << "U. M. Ascher and L. R. Petzold\n"
                << "p. 106\n"
                << "gamma = (2+-sqrt(2))/2\n"
                << "c = [  gamma   1     ]'\n"
                << "A = [  gamma   0     ]\n"
                << "    [ 1-gamma  gamma ]\n"
                << "b = [ 1-gamma  gamma ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    pl->set<double>("gamma",gammaDefault_,
      "The default value is gamma = (2-sqrt(2))/2. "
      "This will produce an L-stable 2nd order method with the stage "
      "times within the timestep.  Other values of gamma will still "
      "produce an L-stable scheme, but will only be 1st order accurate.");

    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) =                              gamma_; A(0,1) = zero;
    A(1,0) = Teuchos::as<Scalar>( one - gamma_ ); A(1,1) = gamma_;

    // Fill b:
    b(0) = Teuchos::as<Scalar>( one - gamma_ ); b(1) = gamma_;

    // Fill c:
    c(0) = gamma_; c(1) = one;

    int order = 1;
    if ( std::abs((gamma_-gammaDefault_)/gamma_) < 1.0e-08 ) order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }

  private:
    Scalar gammaDefault_;
    Scalar gamma_;
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 3 Stage 2nd  order
 *
 *  The tableau (order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}  \gamma  & \gamma      &         &        \\
 *                       1-\gamma & 1-2\gamma   & \gamma  &        \\
 *                       1-2      & 1/2 -\gamma & 0       & \gamma \\ \hline
 *                                & 1/6         & 1/6     & 2/3
 *  \end{array}
 *  \f]
 *  The value is \f$\gamma = 1/ (2 + \sqrt{2})\f$.
 *  This will produce an L-stable 2nd order method with the stage
 *  times within the timestep.
 *
 *  Reference: Implicit-explicit Runge-Kutta schemes and applications to
 *             hyperbolic systems with relaxation
 *             L Pareschi, G Russo
 *             Journal of Scientific computing, 2005 - Springer
 *             Table 5
 */

template<class Scalar>
class StepperSDIRK_3Stage2ndOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_3Stage2ndOrder()
  {
    this->setStepperType("SDIRK 3 Stage 2nd order");
    this->setupTableau();
    this->setupDefault();
  }
#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_3Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {

    this->setStepperType("SDIRK 3 Stage 2nd order");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_3Stage2ndOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {

    this->setStepperType("SDIRK 3 Stage 2nd order");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Implicit-explicit Runge-Kutta schemes and applications to\n"
                << "hyperbolic systems with relaxation\n"
                << "L Pareschi, G Russo\n"
                << "Journal of Scientific computing, 2005 - Springer\n"
                << "Table 5\n"
                << "gamma = 1/(2+sqrt(2))\n"
                << "c = [  gamma   (1-gamma)   1/2  ]'\n"
                << "A = [  gamma      0         0   ]\n"
                << "    [ 1-2gamma   gamma      0   ]\n"
                << "    [ 1/2-gamma   0      gamma  ]\n"
                << "b = [  1/6       1/6       2/3  ]'";
    return Description.str();
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 3;
    const int order = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar gamma = as<Scalar>(one - ( one / ST::squareroot(2*one) ) );

    // Fill A:
    A(0,0) = A(1,1) = A(2,2) = gamma;
    A(0,1) = A(0,2) = A(1,2) = A(2,1) = zero;
    A(1,0) = as<Scalar>(one - 2*gamma);
    A(2,0) = as<Scalar>(  ( one/ (2.*one)) - gamma );

    // Fill b:
    b(0) = b(1) = ( one / (6*one) );
    b(2) = (2*one)/(3*one);

    // Fill c:
    c(0) = gamma;
    c(1) = one - gamma;
    c(2) = one / (2*one);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }

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
 *  The default value is \f$\gamma = (3 + \sqrt{3})/6\f$.
 *
 *  Reference: E. Hairer, S. P. Norsett, and G. Wanner,
 *             Solving Ordinary Differential Equations I:
 *             Nonstiff Problems, 2nd Revised Edition,
 *             Table 7.2, pg 207.
 */
template<class Scalar>
class StepperSDIRK_2Stage3rdOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_2Stage3rdOrder()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const Scalar one = ST::one();
    gammaDefault_ = as<Scalar>((3*one+ST::squareroot(3*one))/(6*one));
    gammaTypeDefault_ = "3rd Order A-stable";

    this->setStepperType("SDIRK 2 Stage 3rd order");
    this->setGammaType(gammaTypeDefault_);
    this->setGamma(gammaDefault_);
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    std::string gammaType = "3rd Order A-stable",
    Scalar gamma = Scalar(0.7886751345948128))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const Scalar one = ST::one();
    gammaDefault_ = as<Scalar>((3*one+ST::squareroot(3*one))/(6*one));
    gammaTypeDefault_ = "3rd Order A-stable";

    this->setStepperType("SDIRK 2 Stage 3rd order");
    this->setGammaType(gammaType);
    this->setGamma(gamma);
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    std::string gammaType = "3rd Order A-stable",
    Scalar gamma = Scalar(0.7886751345948128))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const Scalar one = ST::one();
    gammaDefault_ = as<Scalar>((3*one+ST::squareroot(3*one))/(6*one));
    gammaTypeDefault_ = "3rd Order A-stable";

    this->setStepperType("SDIRK 2 Stage 3rd order");
    this->setGammaType(gammaType);
    this->setGamma(gamma);
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  void setGammaType(std::string gammaType)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !(gammaType == "3rd Order A-stable" or
        gammaType == "2nd Order L-stable" or
        gammaType == "gamma"), std::logic_error,
      "gammaType needs to be '3rd Order A-stable', '2nd Order L-stable' "
      "or 'gamma'.");

    gammaType_ = gammaType;
    this->isInitialized_ = false;
    this->setupTableau();
  }

  std::string getGammaType() { return gammaType_; }

  void setGamma(Scalar gamma)
  {
    if ( gammaType_ == "gamma" ) {
      gamma_ = gamma;
      this->setupTableau();
    }
    this->isInitialized_ = false;
  }

  Scalar getGamma() { return gamma_; }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.2, pg 207\n"
                << "gamma = (3+sqrt(3))/6 -> 3rd order and A-stable\n"
                << "gamma = (2-sqrt(2))/2 -> 2nd order and L-stable\n"
                << "c = [  gamma     1-gamma  ]'\n"
                << "A = [  gamma     0        ]\n"
                << "    [ 1-2*gamma  gamma    ]\n"
                << "b = [ 1/2        1/2      ]'";
    return Description.str();
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);

    pl->set<bool>("Initial Condition Consistency Check", false);
    pl->set<std::string>("Gamma Type", gammaTypeDefault_,
      "Valid values are '3rd Order A-stable' ((3+sqrt(3))/6.), "
      "'2nd Order L-stable' ((2-sqrt(2))/2) and 'gamma' (user defined).  "
      "The default value is '3rd Order A-stable'.");
    pl->set<double>("gamma", gammaDefault_,
      "Equal to (3+sqrt(3))/6 if 'Gamma Type' = '3rd Order A-stable', or "
      "(2-sqrt(2))/2 if 'Gamma Type' = '2nd Order L-stable', or "
      "user-defined gamma value if 'Gamma Type = 'gamma'.  "
      "The default value is gamma = (3+sqrt(3))/6, which matches "
      "the default 'Gamma Type' = '3rd Order A-stable'.");

    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    int order = 0;
    if (gammaType_ == "3rd Order A-stable") {
      order = 3;
      gamma_ = gammaDefault_;
    } else if (gammaType_ == "2nd Order L-stable") {
      order = 2;
      gamma_ = as<Scalar>( (2*one - ST::squareroot(2*one))/(2*one) );
    } else if (gammaType_ == "gamma") {
      order = 2;
    }

    // Fill A:
    A(0,0) =                     gamma_; A(0,1) = zero;
    A(1,0) = as<Scalar>(one - 2*gamma_); A(1,1) = gamma_;

    // Fill b:
    b(0) = as<Scalar>( one/(2*one) ); b(1) = as<Scalar>( one/(2*one) );

    // Fill c:
    c(0) = gamma_; c(1) = as<Scalar>( one - gamma_ );

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,2,3));
  }

  private:
    std::string gammaTypeDefault_;
    std::string gammaType_;
    Scalar      gammaDefault_;
    Scalar      gamma_;
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
class StepperEDIRK_2Stage3rdOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperEDIRK_2Stage3rdOrder()
  {
    this->setStepperType("EDIRK 2 Stage 3rd order");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperEDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("EDIRK 2 Stage 3rd order");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperEDIRK_2Stage3rdOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("EDIRK 2 Stage 3rd order");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Hammer & Hollingsworth method\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.1, pg 205\n"
                << "c = [  0   2/3 ]'\n"
                << "A = [  0    0  ]\n"
                << "    [ 1/3  1/3 ]\n"
                << "b = [ 1/4  3/4 ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    // Fill A:
    A(0,0) =                      zero; A(0,1) =                      zero;
    A(1,0) = as<Scalar>( one/(3*one) ); A(1,1) = as<Scalar>( one/(3*one) );

    // Fill b:
    b(0) = as<Scalar>( one/(4*one) ); b(1) = as<Scalar>( 3*one/(4*one) );

    // Fill c:
    c(0) = zero; c(1) = as<Scalar>( 2*one/(3*one) );
    int order = 3;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 1 Stage Theta
 *
 *  The tableau (order = 1 or 2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} \theta  & \theta  \\
 *                             & 1       \end{array}
 *  \f]
 *  Valid values are \f$0 < \theta \leq 1\f$, where \f$\theta\f$ = 0
 *  implies Forward Euler (not avialble with this stepepr as it makes it
 *  explicit), \f$theta\f$ = 1/2 implies implicit midpoint
 *  method (default), and \f$theta\f$ = 1 implies Backward Euler.
 *  For \f$theta\f$ != 1/2, this method is first-order accurate,
 *  and with \f$theta\f$ = 1/2, it is second-order accurate.
 *  This method is A-stable, but becomes L-stable with \f$theta\f$ = 1.
 *  (A.K.A. Generalized Implicit Midpoint Method.)
 *
 *  Reference: Non-standard finite-difference methods
 *             in dynamical systems, P. Kama,
 *             Dissertation, University of Pretoria, pg. 49.
 */
template<class Scalar>
class StepperDIRK_1StageTheta :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperDIRK_1StageTheta()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("DIRK 1 Stage Theta Method");
    this->setTheta(thetaDefault_);
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperDIRK_1StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    Scalar theta = Scalar(0.5))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("DIRK 1 Stage Theta Method");
    this->setTheta(theta);
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperDIRK_1StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    Scalar theta = Scalar(0.5))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("DIRK 1 Stage Theta Method");
    this->setTheta(theta);
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  void setTheta(Scalar theta)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      theta == Teuchos::ScalarTraits<Scalar>::zero(), std::logic_error,
      "'theta' can not be zero, as it makes this stepper explicit. \n"
      "Try using the 'RK Forward Euler' stepper.\n");
    theta_ = theta;
    this->setupTableau();
    this->isInitialized_ = false;
  }

  Scalar getTheta() { return theta_; }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Non-standard finite-difference methods\n"
                << "in dynamical systems, P. Kama,\n"
                << "Dissertation, University of Pretoria, pg. 49.\n"
                << "Comment:  Generalized Implicit Midpoint Method\n"
                << "c = [ theta ]'\n"
                << "A = [ theta ]\n"
                << "b = [  1  ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);

    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    pl->set<double>("theta",thetaDefault_,
      "Valid values are 0 <= theta <= 1, where theta = 0 "
      "implies Forward Euler, theta = 1/2 implies implicit midpoint "
      "method (default), and theta = 1 implies Backward Euler. "
      "For theta != 1/2, this method is first-order accurate, "
      "and with theta = 1/2, it is second-order accurate.  "
      "This method is A-stable, but becomes L-stable with theta=1.");

    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    A(0,0) = theta_;
    b(0) = ST::one();
    c(0) = theta_;

    int order = 1;
    if ( std::abs((theta_-thetaDefault_)/theta_) < 1.0e-08 ) order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,1,2));
  }

  private:
    Scalar thetaDefault_;
    Scalar theta_;
};


// ----------------------------------------------------------------------------
/** \brief EDIRK 2 Stage Theta Method
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0       &        \\
 *                       1  & 1-\theta & \theta \\ \hline
 *                          & 1-\theta & \theta  \end{array}
 *  \f]
 *  Valid values are \f$0 < \theta <= 1\f$, where \f$\theta\f$ = 0
 *  implies Forward Euler (not avialble with this stepepr as it makes it
 *  explicit), \f$\theta\f$ = 1/2 implies trapezoidal
 *  method (default), and \f$\theta\f$ = 1 implies Backward Euler.
 *  For \f$\theta\f$ != 1/2, this method is first-order accurate,
 *  and with \f$\theta\f$ = 1/2, it is second-order accurate.
 *  This method is A-stable, but becomes L-stable with \f$\theta\f$=1.
 *
 *  Reference: Computer Methods for ODEs and DAEs,
 *             U. M. Ascher and L. R. Petzold, p. 113.
 */
template<class Scalar>
class StepperEDIRK_2StageTheta :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperEDIRK_2StageTheta()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("EDIRK 2 Stage Theta Method");
    this->setTheta(thetaDefault_);
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperEDIRK_2StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    Scalar theta = Scalar(0.5))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("EDIRK 2 Stage Theta Method");
    this->setTheta(theta);
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperEDIRK_2StageTheta(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    Scalar theta = Scalar(0.5))
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    thetaDefault_ = ST::one()/(2*ST::one());

    this->setStepperType("EDIRK 2 Stage Theta Method");
    this->setTheta(theta);
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  void setTheta(Scalar theta)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      theta == Teuchos::ScalarTraits<Scalar>::zero(), std::logic_error,
      "'theta' can not be zero, as it makes this stepper explicit. \n"
      "Try using the 'RK Forward Euler' stepper.\n");
    theta_ = theta;
    this->isInitialized_ = false;
    this->setupTableau();
  }

  Scalar getTheta() { return theta_; }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Computer Methods for ODEs and DAEs\n"
                << "U. M. Ascher and L. R. Petzold\n"
                << "p. 113\n"
                << "c = [  0       1     ]'\n"
                << "A = [  0       0     ]\n"
                << "    [ 1-theta  theta ]\n"
                << "b = [ 1-theta  theta ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);

    pl->set<bool>("Initial Condition Consistency Check", false);
    pl->set<double>("theta",thetaDefault_,
      "Valid values are 0 < theta <= 1, where theta = 0 "
      "implies Forward Euler, theta = 1/2 implies trapezoidal "
      "method (default), and theta = 1 implies Backward Euler. "
      "For theta != 1/2, this method is first-order accurate, "
      "and with theta = 1/2, it is second-order accurate.  "
      "This method is A-stable, but becomes L-stable with theta=1.");

    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =                               zero; A(0,1) =  zero;
    A(1,0) = Teuchos::as<Scalar>( one - theta_ ); A(1,1) = theta_;

    // Fill b:
    b(0) = Teuchos::as<Scalar>( one - theta_ );
    b(1) = theta_;

    // Fill c:
    c(0) = zero;
    c(1) = one;

    int order = 1;
    if ( std::abs((theta_-thetaDefault_)/theta_) < 1.0e-08 ) order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,1,2));
  }

  private:
    Scalar thetaDefault_;
    Scalar theta_;
};


// ----------------------------------------------------------------------------
/** \brief RK Trapezoidal Rule (A.K.A. RK Crank-Nicolson)
 *
 *  The tableau (order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  &  0  &     \\
 *                       1  & 1/2 & 1/2 \\ \hline
 *                          & 1/2 & 1/2  \end{array}
 *  \f]
 *  It is second-order accurate and A-stable.
 *
 *  Reference: Computer Methods for ODEs and DAEs,
 *             U. M. Ascher and L. R. Petzold, p. 113.
 */
template<class Scalar>
class StepperEDIRK_TrapezoidalRule :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperEDIRK_TrapezoidalRule()
  {
    this->setStepperType("RK Trapezoidal Rule");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperEDIRK_TrapezoidalRule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("RK Trapezoidal Rule");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperEDIRK_TrapezoidalRule(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Trapezoidal Rule");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Also known as Crank-Nicolson Method.\n"
                << "c = [  0   1   ]'\n"
                << "A = [  0   0   ]\n"
                << "    [ 1/2  1/2 ]\n"
                << "b = [ 1/2  1/2 ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    const Scalar one = ST::one();
    const Scalar zero = ST::zero();
    const Scalar onehalf = ST::one()/(2*ST::one());

    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    // Fill A:
    A(0,0) =    zero; A(0,1) =    zero;
    A(1,0) = onehalf; A(1,1) = onehalf;

    // Fill b:
    b(0) = onehalf;
    b(1) = onehalf;

    // Fill c:
    c(0) = zero;
    c(1) = one;

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK Implicit Midpoint
 *
 *  The tableau (order = 1 or 2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c} 1/2  & 1/2  \\ \hline
 *                          & 1    \end{array}
 *  \f]
 *  Implicit midpoint method is second-order accurate, and is A-stable.
 *
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner,
 *             Table 5.2, pg 72.
 *
 *             Solving Ordinary Differential Equations I:
 *             Nonstiff Problems, 2nd Revised Edition,
 *             E. Hairer, S. P. Norsett, and G. Wanner,
 *             Table 7.1, pg 205,
 */
template<class Scalar>
class StepperSDIRK_ImplicitMidpoint :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_ImplicitMidpoint()
  {
    this->setStepperType("RK Implicit Midpoint");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_ImplicitMidpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("RK Implicit Midpoint");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_ImplicitMidpoint(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Implicit Midpoint");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.2, pg 72\n"
                << "Solving Ordinary Differential Equations I:\n"
                << "Nonstiff Problems, 2nd Revised Edition\n"
                << "E. Hairer, S. P. Norsett, and G. Wanner\n"
                << "Table 7.1, pg 205\n"
                << "c = [ 1/2 ]'\n"
                << "A = [ 1/2 ]\n"
                << "b = [  1  ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    int NumStages = 1;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar onehalf = ST::one()/(2*ST::one());
    const Scalar one = ST::one();

    // Fill A:
    A(0,0) = onehalf;

    // Fill b:
    b(0) = one;

    // Fill c:
    c(0) = onehalf;

    int order = 2;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Diagonally-Implicit RK Butcher Tableau
 *
 *  The tableau (stage=2, order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  1/4  & 1/4  & \\
 *                         3/4  & 1/2  & 1/4 \\   \hline
 *                              & 1/2  & 1/2  \end{array}
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 */
template<class Scalar>
class StepperSDIRK_SSPDIRK22 :
  virtual public StepperDIRK<Scalar>
{
  public:
  StepperSDIRK_SSPDIRK22()
  {
    this->setStepperType("SSPDIRK22");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_SSPDIRK22(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SSPDIRK22");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_SSPDIRK22(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SSPDIRK22");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
      << "Strong Stability Preserving Diagonally-Implicit RK (stage=2, order=2)\n"
      << "SSP-Coef = 4\n"
      << "c =     [ 1/4   3/4 ]'\n"
      << "A =     [ 1/4       ]\n"
      << "        [ 1/2   1/4 ]\n"
      << "b     = [ 1/2   1/2 ]\n" << std::endl;
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 2;
    const int order     = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one       = ST::one();
    const Scalar zero      = ST::zero();
    const Scalar onehalf   = one/(2*one);
    const Scalar onefourth = one/(4*one);

    // Fill A:
    A(0,0) = A(1,1) = onefourth;
    A(0,1) = zero;
    A(1,0) = onehalf;

    // Fill b:
    b(0) = b(1) = onehalf;

    // Fill c:
    c(0) = A(0,0);
    c(1) = A(1,0) + A(1,1);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Diagonally-Implicit RK Butcher Tableau
 *
 *  The tableau (stage=3, order=2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  1/6  & 1/6  &           \\
 *                         1/2  & 1/3  & 1/6 &     \\
 *                         5/6  & 1/3  & 1/3 & 1/3 \\   \hline
 *                              & 1/3  & 1/3 & 1/3  \end{array}
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 */
template<class Scalar>
class StepperSDIRK_SSPDIRK32 :
  virtual public StepperDIRK<Scalar>
{
  public:
  StepperSDIRK_SSPDIRK32()
  {
    this->setStepperType("SSPDIRK32");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_SSPDIRK32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SSPDIRK32");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_SSPDIRK32(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SSPDIRK32");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "Strong Stability Preserving Diagonally-Implicit RK (stage=3, order=2)\n"
                << "SSP-Coef = 6\n"
                << "c = [ 1/6   1/2   5/6 ]'\n"
                << "A = [ 1/6             ]\n"
                << "    [ 1/3   1/6       ]\n"
                << "    [ 1/3   1/3   1/6 ]\n"
                << "b = [ 1/3   1/3   1/3 ]\n" << std::endl;
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 3;
    const int order     = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one      = ST::one();
    const Scalar zero     = ST::zero();
    const Scalar onethird = one/(3*one);
    const Scalar onesixth = one/(6*one);

    // Fill A:
    A(0,0) = A(1,1) = A(2,2) = onesixth;
    A(1,0) = A(2,0) = A(2,1) = onethird;
    A(0,1) = A(0,2) = A(1,2) = zero;

    // Fill b:
    b(0) = b(1) = b(2) = onethird;

    // Fill c:
    c(0) = A(0,0);
    c(1) = A(1,0) + A(1,1);
    c(2) = A(2,0) + A(2,1) + A(2,2);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Diagonally-Implicit RK Butcher Tableau
 *
 *  The tableau (stage=2, order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  1/(3+\sqrt{3})    & 1/(3+\sqrt{3}) & \\
 *                         (1/6)(3+\sqrt{3}) & 1/\sqrt{3}     & 1/(3+\sqrt{3}) \\   \hline
 *                                           & 1/2            & 1/2 \end{array}
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 */
template<class Scalar>
class StepperSDIRK_SSPDIRK23 :
  virtual public StepperDIRK<Scalar>
{
  public:
  StepperSDIRK_SSPDIRK23()
  {
    this->setStepperType("SSPDIRK23");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_SSPDIRK23(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SSPDIRK23");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_SSPDIRK23(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SSPDIRK23");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
      << "Strong Stability Preserving Diagonally-Implicit RK (stage=2, order=3)\n"
      << "SSP-Coef = 1 + sqrt( 3 )\n"
      << "c =     [ 1/(3 + sqrt( 3 ))  (1/6)(3 + sqrt( 3 )) ] '\n"
      << "A =     [ 1/(3 + sqrt( 3 ))                       ] \n"
      << "        [ 1/sqrt( 3 )        1/(3 + sqrt( 3 ))    ] \n"
      << "b     = [ 1/2                   1/2               ] \n" << std::endl;
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 2;
    const int order     = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one       = ST::one();
    const Scalar zero      = ST::zero();
    const Scalar onehalf   = one/(2*one);
    const Scalar rootthree = ST::squareroot(3*one);

    // Fill A:
    A(0,0) = A(1,1) = one/(3*one + rootthree);
    A(1,0) = one/rootthree;
    A(0,1) = zero;

    // Fill b:
    b(0) = b(1) = onehalf;

    // Fill c:
    c(0) = A(0,0);
    c(1) = A(1,0) + A(1,1);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief Strong Stability Preserving Diagonally-Implicit RK Butcher Tableau
 *
 *  The tableau (stage=3, order=3) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T \\
 *      & \hat{b}^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  1/6  & 1/6  & \\
 *                         1/2  & 1/3  & 1/6 & \\
 *                         5/6  & 1/3  & 1/3 & 1/3 \\   \hline
 *                              & 1/3  & 1/3 & 1/3  \end{array}
 *  \f]
 *  Reference:  Gottlieb, S., Ketcheson, D.I., Shu, C.-W.
 *              Strong Stability Preserving Runge–Kutta and Multistep Time Discretizations.
 *              World Scientific Press, London (2011)
 */
template<class Scalar>
class StepperSDIRK_SSPDIRK33 :
  virtual public StepperDIRK<Scalar>
{
  public:
  StepperSDIRK_SSPDIRK33()
  {
    this->setStepperType("SSPDIRK33");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_SSPDIRK33(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SSPDIRK33");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_SSPDIRK33(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SSPDIRK33");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
      << "Strong Stability Preserving Diagonally-Implicit RK (stage=3, order=3)\n"
      << "SSP-Coef = 2 + 2 sqrt(2)\n"
      << "c =     [ 1/( 4 + 2 sqrt(2)      1/2            (1/4)(2 + sqrt(2) ] '\n"
      << "A =     [ 1/( 4 + 2 sqrt(2)                                       ] \n"
      << "        [ 1/(2 sqrt(2)       1/( 4 + 2 sqrt(2)                    ] \n"
      << "        [ 1/(2 sqrt(2)        1/(2 sqrt(2)      1/( 4 + 2 sqrt(2) ] \n"
      << "b     = [ 1/3                    1/3                1/3           ] \n"
      << std::endl;
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {

    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    const int NumStages = 3;
    const int order     = 3;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);

    const Scalar one      = ST::one();
    const Scalar zero     = ST::zero();
    const Scalar onethird = one/(3*one);
    const Scalar rootwo   = ST::squareroot(2*one);

    // Fill A:
    A(0,0) = A(1,1) = A(2,2) = one / (4*one + 2*rootwo);
    A(1,0) = A(2,0) = A(2,1) = one / (2*rootwo);
    A(0,1) = A(0,2) = A(1,2) = zero;

    // Fill b:
    b(0) = b(1) = b(2) = onethird;

    // Fill c:
    c(0) = A(0,0);
    c(1) = A(1,0) + A(1,1);
    c(2) = A(2,0) + A(2,1) + A(2,2);

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief RK Implicit 1 Stage 1st order Radau IA
 *
 *  The tableau (order = 1) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|c}  0  & 1  \\ \hline
 *                         & 1  \end{array}
 *  \f]
 *  and is A-stable.
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner,
 *             Table 5.3, pg 73.
 */
template<class Scalar>
class StepperDIRK_1Stage1stOrderRadauIA :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperDIRK_1Stage1stOrderRadauIA()
  {
    this->setStepperType("RK Implicit 1 Stage 1st order Radau IA");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperDIRK_1Stage1stOrderRadauIA(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("RK Implicit 1 Stage 1st order Radau IA");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperDIRK_1Stage1stOrderRadauIA(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Implicit 1 Stage 1st order Radau IA");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.3, pg 73\n"
                << "c = [ 0 ]'\n"
                << "A = [ 1 ]\n"
                << "b = [ 1 ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
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

    auto emptyBStar = Teuchos::SerialDenseVector<int,Scalar>();
    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,emptyBStar,false));
  }
};


// ----------------------------------------------------------------------------
/** \brief RK Implicit 2 Stage 2nd order Lobatto IIIB
 *
 *  The tableau (order = 2) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cc}  0  & 1/2 &  0  \\
 *                       1  & 1/2 &  0  \\ \hline
 *                          & 1/2 & 1/2  \end{array}
 *  \f]
 *  It is second-order accurate and A-stable.
 *
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner,
 *             Table 5.9, pg 76.
 */
template<class Scalar>
class StepperDIRK_2Stage2ndOrderLobattoIIIB :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperDIRK_2Stage2ndOrderLobattoIIIB()
  {
    this->setStepperType("RK Implicit 2 Stage 2nd order Lobatto IIIB");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperDIRK_2Stage2ndOrderLobattoIIIB(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("RK Implicit 2 Stage 2nd order Lobatto IIIB");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperDIRK_2Stage2ndOrderLobattoIIIB(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("RK Implicit 2 Stage 2nd order Lobatto IIIB");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "Table 5.9, pg 76\n"
                << "c = [  0    1   ]'\n"
                << "A = [ 1/2   0   ]\n"
                << "    [ 1/2   0   ]\n"
                << "b = [ 1/2  1/2  ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 2;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();

    // Fill A:
    A(0,0) = as<Scalar>( one/(2*one) );
    A(0,1) = zero;
    A(1,0) = as<Scalar>( one/(2*one) );
    A(1,1) = zero;

    // Fill b:
    b(0) = as<Scalar>( one/(2*one) );
    b(1) = as<Scalar>( one/(2*one) );

    // Fill c:
    c(0) = zero;
    c(1) = one;
    int order = 2;

    auto emptyBStar = Teuchos::SerialDenseVector<int,Scalar>();
    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,emptyBStar,false));
  }

};


// ----------------------------------------------------------------------------
/** \brief SDIRK 5 Stage 4th order
 *
 *  The tableau (order = 4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}
 *    1/4   & 1/4       &            &         &         &     \\
 *    3/4   & 1/2       & 1/4        &         &         &     \\
 *    11/20 & 17/50     & -1/25      & 1/4     &         &     \\
 *    1/2   & 371/1360  & -137/2720  & 15/544  & 1/4     &     \\
 *     1    & 25/24     & -49/48     & 125/16  & -85/12  & 1/4 \\ \hline
 *          & 25/24     & -49/48     & 125/16  & -85/12  & 1/4 \end{array}
 *  \f]
 *  and is L-stable.
 *
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner, pg100.
 */
template<class Scalar>
class StepperSDIRK_5Stage4thOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_5Stage4thOrder()
  {
    this->setStepperType("SDIRK 5 Stage 4th order");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_5Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SDIRK 5 Stage 4th order");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_5Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SDIRK 5 Stage 4th order");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
      << "L-stable\n"
      << "Solving Ordinary Differential Equations II:\n"
      << "Stiff and Differential-Algebraic Problems,\n"
      << "2nd Revised Edition\n"
      << "E. Hairer and G. Wanner\n"
      << "pg100 \n"
      << "c     = [ 1/4       3/4        11/20   1/2     1   ]'\n"
      << "A     = [ 1/4                                      ]\n"
      << "        [ 1/2       1/4                            ]\n"
      << "        [ 17/50     -1/25      1/4                 ]\n"
      << "        [ 371/1360  -137/2720  15/544  1/4         ]\n"
      << "        [ 25/24     -49/48     125/16  -85/12  1/4 ]\n"
      << "b     = [ 25/24     -49/48     125/16  -85/12  1/4 ]'";
      // << "b     = [ 59/48     -17/96     225/32  -85/12  0   ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::as;
    int NumStages = 5;
    Teuchos::SerialDenseMatrix<int,Scalar> A(NumStages,NumStages);
    Teuchos::SerialDenseVector<int,Scalar> b(NumStages);
    Teuchos::SerialDenseVector<int,Scalar> c(NumStages);
    const Scalar zero = ST::zero();
    const Scalar one = ST::one();
    const Scalar onequarter = as<Scalar>( one/(4*one) );

    // Fill A:
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

    // Fill b:
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

    // Fill c:
    c(0) = onequarter;
    c(1) = as<Scalar>( 3*one/(4*one) );
    c(2) = as<Scalar>( 11*one/(20*one) );
    c(3) = as<Scalar>( one/(2*one) );
    c(4) = one;

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 3 Stage 4th order
 *
 *  The tableau (order = 4) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccc}
 *    \gamma   &     \gamma &           &        \\
 *    1/2      & 1/2-\gamma &    \gamma &        \\
 *    1-\gamma &    2\gamma & 1-4\gamma & \gamma \\ \hline
 *             &     \delta & 1-2\delta & \delta \end{array}
 *  \f]
 *  where \f$\gamma = (1/\sqrt{3})\cos(\pi/18)+1/2\f$ and
 *  \f$\delta = 1/(6(2\gamma-1)^2)\f$, and is A-stable.
 *
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner, p. 100.
 */
template<class Scalar>
class StepperSDIRK_3Stage4thOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_3Stage4thOrder()
  {
    this->setStepperType("SDIRK 3 Stage 4th order");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_3Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SDIRK 3 Stage 4th order");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_3Stage4thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SDIRK 3 Stage 4th order");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "A-stable\n"
                << "Solving Ordinary Differential Equations II:\n"
                << "Stiff and Differential-Algebraic Problems,\n"
                << "2nd Revised Edition\n"
                << "E. Hairer and G. Wanner\n"
                << "p. 100 \n"
                << "gamma = (1/sqrt(3))*cos(pi/18)+1/2\n"
                << "delta = 1/(6*(2*gamma-1)^2)\n"
                << "c = [ gamma      1/2        1-gamma ]'\n"
                << "A = [ gamma                         ]\n"
                << "    [ 1/2-gamma  gamma              ]\n"
                << "    [ 2*gamma    1-4*gamma  gamma   ]\n"
                << "b = [ delta      1-2*delta  delta   ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
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

    // Fill A:
    A(0,0) = gamma;
    A(0,1) = zero;
    A(0,2) = zero;

    A(1,0) = as<Scalar>( one/(2*one) - gamma );
    A(1,1) = gamma;
    A(1,2) = zero;

    A(2,0) = as<Scalar>( 2*gamma );
    A(2,1) = as<Scalar>( one - 4*gamma );
    A(2,2) = gamma;

    // Fill b:
    b(0) = delta;
    b(1) = as<Scalar>( one-2*delta );
    b(2) = delta;

    // Fill c:
    c(0) = gamma;
    c(1) = as<Scalar>( one/(2*one) );
    c(2) = as<Scalar>( one - gamma );

    int order = 4;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 5 Stage 5th order
 *
 *  The tableau (order = 5) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A \\ \hline
 *      & b^T
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|ccccc}
 *    (6-\sqrt{6})/10   & (6-\sqrt{6})/10              & 0                           & 0                       & 0                     & 0                \\
 *    (6+9\sqrt{6})/35  & (-6+5\sqrt{6})/14            & (6-\sqrt{6})/10             & 0                       & 0                     & 0                \\
 *         1            & (888+607\sqrt{6})/2850       & (126-161\sqrt{6})/1425      & (6-\sqrt{6})/10         & 0                     & 0                \\
 *    (4-\sqrt{6})/10   & (3153-3082\sqrt{6})/14250    & (3213+1148\sqrt{6})/28500   & (-267+88\sqrt{6})/500   & (6-\sqrt{6})/10       & 0                \\
 *    (4+\sqrt{6})/10   & (-32583+14638\sqrt{6})/71250 & (-17199+364\sqrt{6})/142500 & (1329-544\sqrt{6})/2500 & (-96+131\sqrt{6})/625 & (6-\sqrt{6})/10  \\ \hline
 *                      &       0                      &       0                     &           1/9           & (16-\sqrt{6})/36      & (16+\sqrt{6})/36
 *  \end{array}
 *  \f]
 *
 *  Reference: Solving Ordinary Differential Equations II:
 *             Stiff and Differential-Algebraic Problems,
 *             2nd Revised Edition, E. Hairer and G. Wanner, pg101.
 */
template<class Scalar>
class StepperSDIRK_5Stage5thOrder :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_5Stage5thOrder()
  {
    this->setStepperType("SDIRK 5 Stage 5th order");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_5Stage5thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SDIRK 5 Stage 5th order");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_5Stage5thOrder(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SDIRK 5 Stage 5th order");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
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
      << "    [ (16+sqrt(6))/36 ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
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

    // Fill A:
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

    // Fill b:
    b(0) = zero;
    b(1) = zero;
    b(2) = as<Scalar>( one/(9*one) );
    b(3) = as<Scalar>( (16*one-sqrt6)/(36*one) );
    b(4) = as<Scalar>( (16*one+sqrt6)/(36*one) );

    // Fill c:
    c(0) = gamma;
    c(1) = as<Scalar>( (6*one+9*one*sqrt6)/(35*one) );
    c(2) = one;
    c(3) = as<Scalar>( (4*one-sqrt6)/(10*one) );
    c(4) = as<Scalar>( (4*one+sqrt6)/(10*one) );

    int order = 5;

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order));
  }
};


// ----------------------------------------------------------------------------
/** \brief SDIRK 2(1) pair
 *
 *  The tableau (order=2(1)) is
 *  \f[
 *  \begin{array}{c|c}
 *    c & A   \\ \hline
 *      & b^T \\
 *      & b^{*T}
 *  \end{array}
 *  \;\;\;\;\mbox{ where }\;\;\;\;
 *  \begin{array}{c|cccc}  0 & 0   & \\
 *                         1 & -1  & 1 \\ \hline
 *                           & 1/2 & 1/2 \\
 *                           & 1   & 0 \end{array}
 *  \f]
 */
template<class Scalar>
class StepperSDIRK_21Pair :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperSDIRK_21Pair()
  {
    this->setStepperType("SDIRK 2(1) Pair");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperSDIRK_21Pair(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess)
  {
    this->setStepperType("SDIRK 2(1) Pair");
    this->setupTableau();
    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperSDIRK_21Pair(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction)
  {
    this->setStepperType("SDIRK 2(1) Pair");
    this->setupTableau();
    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::ostringstream Description;
    Description << this->getStepperType() << "\n"
                << "c =     [  1  0   ]'\n"
                << "A =     [  1      ]\n"
                << "        [ -1  1   ]\n"
                << "b     = [ 1/2 1/2 ]'\n"
                << "bstar = [  1  0   ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());
    return pl;
  }

protected:

  void setupTableau()
  {
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
    A(0,0) =  one; A(0,1) = zero;
    A(1,0) = -one; A(1,1) =  one;

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

    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,order,order,bstar));
  }
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
 *  \begin{array}{c|cc} \gamma  &   \gamma &        \\
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
class StepperDIRK_General :
  virtual public StepperDIRK<Scalar>
{
public:
  /** \brief Default constructor.
   *
   * Requires subsequent setModel() and initialize()
   * calls before calling takestep().
  */
  StepperDIRK_General()
  {
    this->setStepperType("General DIRK");
    this->setupTableau();
    this->setupDefault();
  }

#ifndef TEMPUS_HIDE_DEPRECATED_CODE
  StepperDIRK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<StepperRKObserverComposite<Scalar> >& obs,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    const int order,
    const int orderMin,
    const int orderMax,
    const Teuchos::SerialDenseVector<int,Scalar>& bstar)
  {
    this->setStepperType("General DIRK");
    this->setTableau(A,b,c,order,orderMin,orderMax,bstar);

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->tableau_->isImplicit() != true, std::logic_error,
      "Error - General DIRK did not receive a DIRK Butcher Tableau!\n");

    this->setup(appModel, obs, solver, useFSAL, ICConsistency,
                ICConsistencyCheck, useEmbedded, zeroInitialGuess);
  }
#endif
  StepperDIRK_General(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& appModel,
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
    bool useFSAL,
    std::string ICConsistency,
    bool ICConsistencyCheck,
    bool useEmbedded,
    bool zeroInitialGuess,
    const Teuchos::RCP<StepperRKAppAction<Scalar> >& stepperRKAppAction,
    const Teuchos::SerialDenseMatrix<int,Scalar>& A,
    const Teuchos::SerialDenseVector<int,Scalar>& b,
    const Teuchos::SerialDenseVector<int,Scalar>& c,
    const int order,
    const int orderMin,
    const int orderMax,
    const Teuchos::SerialDenseVector<int,Scalar>& bstar)
  {
    this->setStepperType("General DIRK");
    this->setTableau(A,b,c,order,orderMin,orderMax,bstar);

    TEUCHOS_TEST_FOR_EXCEPTION(
      this->tableau_->isImplicit() != true, std::logic_error,
      "Error - General DIRK did not receive a DIRK Butcher Tableau!\n");

    this->setup(appModel, solver, useFSAL, ICConsistency, ICConsistencyCheck,
                useEmbedded, zeroInitialGuess, stepperRKAppAction);
  }

  std::string getDescription() const
  {
    std::stringstream Description;
    Description << this->getStepperType() << "\n"
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
      << "  gamma = (2-sqrt(2))/2\n"
      << "  c = [  gamma   1     ]'\n"
      << "  A = [  gamma   0     ]\n"
      << "      [ 1-gamma  gamma ]\n"
      << "  b = [ 1-gamma  gamma ]'";
    return Description.str();
  }

  virtual bool getICConsistencyCheckDefault() const { return false; }

  void setupTableau()
  {
    if (this->tableau_ == Teuchos::null) {
      // Set tableau to the default if null, otherwise keep current tableau.
      auto stepper = Teuchos::rcp(new StepperSDIRK_2Stage2ndOrder<Scalar>());
      auto t = stepper->getTableau();
      this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
                                 this->getStepperType(),
                                 t->A(),t->b(),t->c(),
                                 t->order(),t->orderMin(),t->orderMax(),
                                 t->bstar()));
      this->isInitialized_ = false;
    }
  }

  void setTableau(const Teuchos::SerialDenseMatrix<int,Scalar>& A,
                  const Teuchos::SerialDenseVector<int,Scalar>& b,
                  const Teuchos::SerialDenseVector<int,Scalar>& c,
                  const int order,
                  const int orderMin,
                  const int orderMax,
                  const Teuchos::SerialDenseVector<int,Scalar>&
                    bstar = Teuchos::SerialDenseVector<int,Scalar>())
  {
    this->tableau_ = Teuchos::rcp(new RKButcherTableau<Scalar>(
      this->getStepperType(),A,b,c,order,orderMin,orderMax,bstar));
    this->isInitialized_ = false;
  }

  Teuchos::RCP<const Teuchos::ParameterList>
  getValidParameters() const
  {
    Teuchos::RCP<Teuchos::ParameterList> pl = Teuchos::parameterList();
    this->getValidParametersBasicDIRK(pl);
    pl->set<bool>("Initial Condition Consistency Check",
                  this->getICConsistencyCheckDefault());

    // Tableau ParameterList
    Teuchos::RCP<Teuchos::ParameterList> tableauPL = Teuchos::parameterList();
    tableauPL->set<std::string>("A",
     "0.2928932188134524 0.0; 0.7071067811865476 0.2928932188134524");
    tableauPL->set<std::string>("b",
     "0.7071067811865476 0.2928932188134524");
    tableauPL->set<std::string>("c", "0.2928932188134524 1.0");
    tableauPL->set<int>("order", 2);
    tableauPL->set<std::string>("bstar", "");
    pl->set("Tableau", *tableauPL);

    return pl;
  }
};


} // namespace Tempus


#endif // Tempus_StepperRKButcherTableau_hpp
