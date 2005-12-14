#ifndef ARPACK_OPERATORS_HPP
#define ARPACK_OPERATORS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

//! Operators from the ARPACK examples
/*! 
 * These are simple, single processor examples of user-defined
 * Anasazi::Operator-derived classes. The class is templated with ScalarType;
 * possible choices are, for example, "float", "double", 
 * "complex<float>", and "complex<double>"
 *
 * These operators implement those from the ARPACK examples, specifically 
 * covering the example directories:
 * - COMPLEX: nonhermitian complex problems
 * - SYM:     symmetric real problems
 * - NONSYM:  nonsymmetric real problems
 * 
 * Where these problems overlap (e.g., COMPLEX/zndrv1 and NONSYM/dndrv1), 
 * we rely on templating to provide the coverage, instead of implementing
 * the matrices for each data type. Through this mechanism, we are able to 
 * represent the 32 ARPACK examples contained in SYM, NONSYM, and COMPLEX
 * with only 12 sets of operators. 
 * 
 * The following table contains the correspondance between ARPACK examples
 * and the Anasazi operators herein:
 *
 * Anasazi Ex.  | ARPACK Example
 * ---------------------------------------------------------------------
 * ARPACKEx1    | COMPLEX/[zc]ndrv1, NONSYM/[ds]ndrv1
 * ARPACKEx2    | COMPLEX/[zc]ndrv2, NONSYM/[ds]ndrv2
 * ARPACKEx3    | COMPLEX/[zc]ndrv3, NONSYM/[ds]ndrv3
 * ARPACKEx4    | COMPLEX/[zc]ndrv4, NONSYM/[ds]ndrv4
 * ARPACKEx5    | NONSYM/[ds]ndrv5
 * ARPACKEx6    | NONSYM/[ds]ndrv6
 * ARPACKEx7    | SYM/[ds]ndrv1
 * ARPACKEx8    | SYM/[ds]ndrv2
 * ARPACKEx9    | SYM/[ds]ndrv3
 * ARPACKEx10   | SYM/[ds]ndrv4
 * ARPACKEx11   | SYM/[ds]ndrv5
 * ARPACKEx12   | SYM/[ds]ndrv6
 *
 * Because the classes above are templated according to scalar type, 
 * they provide complex versions of the real-only examples from ARPACK.
 * For example, ARPACKEx7 can be templated on real or complex scalar types,
 * although the result should be the same.
 * 
 * The following examples are each contained in a templated struct, with
 * typedefs inside the struct referring to the appropriate classes for
 * the specified operator. For example:
 *   ARPACKEx3<float>::OP  
 * references the operator for the ARPACKEx3 operator, which takes the form 
 * OP = inv[M]*A, templated on the "float" scalar type. Similarly, 
 *   ARPACKEx3<float>::A   and   ARPACKEx3<float>::M 
 * reference the oeprators for the A and M operators, corresponding to 
 * the finite element discretization of the 1-dimensional 
 * convection-diffusion operator (d^2u/dx^2) + rho*(du/dx) on the 
 * interval [0,1] with zero boundary condition using piecewise linear 
 * elements (as specified in COMPLEX/[zc]ndrv3.f and NONSYM/[ds]ndrv3
 * 
 * Documentation for each operator is found with the operator.
 * More documentation is available in
 * packages/anasazi/test/ARPACKExamples/exampledesc
 * 
 * Some code was adapter by the FORTRAN exaples of ARPACK. ARPACK was written 
 * by
 *     Richard Lehoucq
 *     Danny Sorensen
 *     Chao Yang
 *     Dept. of Computational &
 *     Applied Mathematics
 *     Rice University
 *     Houston, Texas
 *
 * See http://www.caam.rice.edu/software/ARPACK/ for more information.
 *
 * \author Chris Baker (SNL/1414,FSU/CSIT) and Heidi Thornquist (SNL/1437)
 *
 * \date Last modified on 12-Dec-05
 */



/*! \class OPA< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the identity
  operator (used in ARPACKEx1, ARPACKEx2, ARPACKEx7, and ARPACKEx8).
*/
template <class ScalarType>
class OPA : public Anasazi::Operator<ScalarType>
{
  
public:
  
  OPA()  {}
  ~OPA() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = Teuchos::ScalarTraits<ScalarType>::zero();
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
      
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
    
    // MyY = ONE*MyX
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
    
    return(Anasazi::Ok);
  }
};



/*! \class OPG< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of 1-D Laplacian
  (like in ARPACKEx7).
*/
template <class ScalarType>
class OPG : public Anasazi::Operator<ScalarType>
{
private:
  void tv(int nx, const ScalarType *x, ScalarType *y) const {
    ScalarType dd = 4.0,
               dl = -Teuchos::ScalarTraits<ScalarType>::one(),
               du = -Teuchos::ScalarTraits<ScalarType>::one();
    
    // Compute the matrix vector multiplication y<---T*x
    // where T is a nx by nx tridiagonal matrix with DD on the 
    // diagonal, DL on the subdiagonal, and DU on the superdiagonal.
    int j;
    j = 0;
    y[j] = dd*x[j] + du*x[j+1];
    for (j=1; j<nx-1; j++) {
      y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
    }
    j = nx-1;
    y[j] = dl*x[j-1] + dd*x[j];
  }
public:
  
  OPG()  {}
  ~OPG() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    int n, nx;
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
      
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
    
    int nvecs = X.GetNumberVecs();
    
    // deduce the size of the operator from the vector length...
    n = X.GetVecLength();
    // ... and the number of interior points in the discretization from that
    nx = Teuchos::ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (nx*nx != n) {
      return(Anasazi::Failed);
    }
    
    // The rest is stolen from the ARPACK codes (see notice above)
    //
    // The matrix used is the 2 dimensional discrete Laplacian on unit
    // square with zero Dirichlet boundary condition.
    //
    // Computes y <--- OP*x, where OP is the nx*nx by nx*nx block 
    // tridiagonal matrix
    //
    //              | T -I          | 
    //              |-I  T -I       |
    //         OP = |   -I  T       |
    //              |        ...  -I|
    //              |           -I T|
    //
    // The subroutine tv is called to computed y<---T*x.
    int p, j, lo;
    for (p=0; p<nvecs; p++) {
      lo = 0;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
    }
    
    for (j=1; j<nx-1; j++) {
      lo = j*nx;
      for (p=0; p<nvecs; p++) {
        tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
        blas.AXPY(nx,-ONE,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
        blas.AXPY(nx,-ONE,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
      }
    }
    
    for (p=0; p<nvecs; p++) {
      lo = (nx-1)*nx;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
    }
    
    return(Anasazi::Ok);
  }
  
};



/*! \class OPE< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  convection-diffusion operator.
*/
template <class ScalarType>
class OPE : public Anasazi::Operator<ScalarType>
{
private:
  ScalarType _rho;
  void tv(int nx, const ScalarType *x, ScalarType *y) const {
    const ScalarType ONE = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType h = ONE / (ScalarType)(nx+1),
               h2 = h*h,
               dd = 4.0 / h2,
               dl = -ONE / h2 - 0.5*_rho / h,
               du = -ONE / h2 + 0.5*_rho / h;
    
    
    // Compute the matrix vector multiplication y<---T*x
    // where T is a nx by nx tridiagonal matrix with DD on the 
    // diagonal, DL on the subdiagonal, and DU on the superdiagonal.
    //
    // When rho*h/2 <= 1, the discrete convection-diffusion operator 
    // has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
    // eigenvalues.
    
    int j;
    j = 0;
    y[j] = dd*x[j] + du*x[j+1];
    for (j=1; j<nx-1; j++) {
      y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
    }
    j = nx-1;
    y[j] = dl*x[j-1] + dd*x[j];
  }
    
public:
  
  OPE(ScalarType rho = Teuchos::ScalarTraits<ScalarType>::zero()) : _rho(rho) {}
  ~OPE() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    int n, nx;
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
      
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
    
    int nvecs = X.GetNumberVecs();
    
    // deduce the size of the operator from the vector length...
    n = X.GetVecLength();
    // ... and the number of interior points in the discretization from that
    nx = Teuchos::ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (nx*nx != n) {
      return(Anasazi::Failed);
    }
    
    // The rest is stolen from the ARPACK codes (see notice above)
    //
    // Computes y <--- OP*x, where OP is the nx*nx by nx*nx block 
    // tridiagonal matrix
    //
    //              | T -I          | 
    //              |-I  T -I       |
    //         OP = |   -I  T       |
    //              |        ...  -I|
    //              |           -I T|
    //
    // derived from the standard central difference discretization 
    // of the 2 dimensional convection-diffusion operator 
    // (Laplacian u) + rho*(du/dx) on a unit square with zero boundary 
    // condition.
    //
    // When rho*h/2 <= 1, the discrete convection-diffusion operator 
    // has real eigenvalues.  When rho*h/2 > 1, it has COMPLEX 
    // eigenvalues.
    //
    // The subroutine TV is called to compute y<---T*x.
    int p, j, lo;
    const ScalarType h2 = ONE / (ScalarType)((nx+1)*(nx+1));
    for (p=0; p<nvecs; p++) {
      lo = 0;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE/h2,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
    }
    
    for (j=1; j<nx-1; j++) {
      lo = j*nx;
      for (p=0; p<nvecs; p++) {
        tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
        blas.AXPY(nx,-ONE/h2,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
        blas.AXPY(nx,-ONE/h2,&(*MyX)[p][lo+nx],1,&(*MyY)[p][lo],1);
      }
    }
    
    for (p=0; p<nvecs; p++) {
      lo = (nx-1)*nx;
      tv(nx,&(*MyX)[p][lo],&(*MyY)[p][lo]);
      blas.AXPY(nx,-ONE/h2,&(*MyX)[p][lo-nx],1,&(*MyY)[p][lo],1);
    }
    
    return(Anasazi::Ok);
  }
};



/*! \class OPF< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  convection-diffusion operator.
  
  The operator applied is:
    OPF = inv[A-sigma*I]
*/
template <class ScalarType>
class OPF : public Anasazi::Operator<ScalarType>
{
private:
  int _n,_nx;
  ScalarType _rho, _sigma;
  std::vector<ScalarType> _dl, _dd, _du, _du2;
  std::vector<int> _ipiv;
  int _ferror;
  
public:
  
  OPF( int n,
       ScalarType rho = Teuchos::ScalarTraits<ScalarType>::zero(),
       ScalarType sigma = Teuchos::ScalarTraits<ScalarType>::zero() ) 
      : _n(n), _rho(rho), _sigma(sigma) {
    
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    const ScalarType ONE = SCT::one();
    const ScalarType TWO = ONE + ONE;
    Teuchos::LAPACK<int,ScalarType> lapack;
    
    _nx = Teuchos::ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (_nx*_nx != n) {
      cout << "Argument 1 to OPF() was not a square number." << endl;
      _n = 100;
      _nx = 10;
    }
    
    /*----------------------------------------------------*
    | Construct C = A - SIGMA*I and factor C using LAPACK |
    | subroutine dgttrf. The matrix A is chosen to be     |
    | the tridiagonal matrix derived from standard        |
    | central difference of the 1-d convection diffusion  |
    | operator u" + rho*u' on the interval [0, 1] with    |
    | zero Dirichlet boundary condition.                  |
    \----------------------------------------------------*/
    ScalarType h,s,s1,s2,s3;
    h = ONE / (ScalarType)(_n+1);
    s = _rho*h / TWO;
    
    s1 = -ONE-s;
    s2 = TWO - SCT::real(sigma);
    s3 = -ONE+s;
    _dl.resize(n-1,s1);
    _dd.resize(n  ,s2);
    _du.resize(n-1,s3);
    _du2.resize(n-2);
    _ipiv.resize(n);
  
    int _ferror;
    lapack.GTTRF(_n,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&_ferror);
    if (_ferror != 0) {
      cout << "Error in GTTRF in OPF()" << endl;
    }
    
  }
  ~OPF() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = Teuchos::ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
    
    // if there were problems with the factorization, quit now
    if (_ferror) {
      return Anasazi::Failed;
    }
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() == Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() == Y.GetVecLength()) return Anasazi::Failed;
    if (X.GetVecLength() == _n) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X = inv[A-SIGMA*I]*X using GTTRS
    int p;
    // set Y = X, as GTTRS operates in situ
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
    // call GTTRS multiple times (it takes multiple RHS, but MyMultiVec doesn't
    // use block storage)
    int ierr;
    for (p=0; p<nvecs; p++) {
      lapack.GTTRS('N',_n,1,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&(*MyY[p][0]),_n,&ierr);
      if (ierr != 0) {
        return Anasazi::Failed;
      }
    }
    
    return(Anasazi::Ok);
  }
};



template <class ScalarType>
struct ARPACKEx1
{
  typedef   OPE<ScalarType>   OP;
  typedef   OPA<ScalarType>   M;
  static const bool isHermitian = false;
};

template <class ScalarType>
struct ARPACKEx2
{
  typedef   OPF<ScalarType>   OP;
  typedef   OPA<ScalarType>   M;
  static const bool isHermitian = false;
};

/*
template <class ScalarType>
struct ARPACKEx3
{
  typedef   OPN<ScalarType>   OP;
  typedef   OPM<ScalarType>   M;
  static const bool isHermitian = false;
};

template <class ScalarType>
struct ARPACKEx4
{
  typedef   OPO<ScalarType>   OP;
  typedef   OPM<ScalarType>   M;
  static const bool isHermitian = false;
};

template <class ScalarType>
struct ARPACKEx5
{
  typedef   OPC<ScalarType>   OP;
  typedef   OPB<ScalarType>   M;
  static const bool isHermitian = false;
};

template <class ScalarType>
struct ARPACKEx6
{
  typedef   OPD<ScalarType>   OP;
  typedef   OPB<ScalarType>   M;
  static const bool isHermitian = false;
};
*/

template <class ScalarType>
struct ARPACKEx7
{
  typedef   OPG<ScalarType>   OP;
  typedef   OPA<ScalarType>   M;
  static const bool isHermitian = true;
};

/*
template <class ScalarType>
struct ARPACKEx8
{
  typedef   OPI<ScalarType>   OP;
  typedef   OPA<ScalarType>   M;
  static const bool isHermitian = true;
};

template <class ScalarType>
struct ARPACKEx9
{
  typedef   OPR<ScalarType>   OP;
  typedef   OPQ<ScalarType>   M;
  static const bool isHermitian = true;
};

template <class ScalarType>
struct ARPACKEx10
{
  typedef   OPS<ScalarType>   OP;
  typedef   OPQ<ScalarType>   M;
  static const bool isHermitian = true;
};

template <class ScalarType>
struct ARPACKEx11
{
  typedef   OPU<ScalarType>   OP;
  typedef   OPP<ScalarType>   M;
  static const bool isHermitian = true;
};

template <class ScalarType>
struct ARPACKEx12
{
  typedef   OPT<ScalarType>   OP;
  typedef   OPQ<ScalarType>   M;
  static const bool isHermitian = true;
};
*/

#endif //ARPACK_OPERATORS_HPP
