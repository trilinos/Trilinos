#ifndef ARPACK_OPERATORS_HPP
#define ARPACK_OPERATORS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;

//! Operators from the ARPACK examples
/*! 
 * These are simple, single processor examples of user-defined
 * Anasazi::Operator-derived classes. The classes are templated with ScalarType;
 * possible choices are, for example, "float", "double", 
 * "complex<float>", and "complex<double>".
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
 * with only 12 sets of operators: 6 symmetric and 6 nonsymmetric.
 * 
 * The following table contains the correspondance between ARPACK examples
 * and the Anasazi operators herein:
 *
 * Anasazi Ex.   | ARPACK example drivers
 * ---------------------------------------------------------------------
 * ARPACK_NDRV1  | COMPLEX/[zc]ndrv1, NONSYM/[ds]ndrv1
 * ARPACK_NDRV2  | COMPLEX/[zc]ndrv2, NONSYM/[ds]ndrv2
 * ARPACK_NDRV3  | COMPLEX/[zc]ndrv3, NONSYM/[ds]ndrv3
 * ARPACK_NDRV4  | COMPLEX/[zc]ndrv4, NONSYM/[ds]ndrv4
 * ARPACK_NDRV5  | NONSYM/[ds]ndrv5
 * ARPACK_NDRV6  | NONSYM/[ds]ndrv6
 * ARPACK_SDRV1  | SYM/[ds]sdrv1
 * ARPACK_SDRV2  | SYM/[ds]sdrv2
 * ARPACK_SDRV3  | SYM/[ds]sdrv3
 * ARPACK_SDRV4  | SYM/[ds]sdrv4
 * ARPACK_SDRV5  | SYM/[ds]sdrv5
 * ARPACK_SDRV6  | SYM/[ds]sdrv6
 *
 * Because the classes above are templated according to scalar type, 
 * they provide complex versions of the real-only examples from ARPACK.
 * For example, ARPACK_SDRV1 can be templated on real or complex scalar types,
 * though the results should be the same.
 * 
 * The following examples are each contained in a templated struct, with
 * typedefs inside the struct referring to the appropriate classes for
 * the specified operator. For example:
 *   ARPACK_NDRV3<float>::OP  
 * references the operator for the ARPACK_NDRV3 operator, which takes the form 
 * OP = inv[M]*A, templated on the "float" scalar type. Similarly, 
 *   ARPACK_NDRV3<float>::A   and   ARPACK_NDRV3<float>::M 
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
 * Some code was adapted from the FORTRAN exaples of ARPACK. ARPACK was written 
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




/******************************************************************************/
/*! \class OPA< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the identity
  operator.
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
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    // MyY = ONE*MyX
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
    
    return(Anasazi::Ok);
  }
};




/******************************************************************************/
/*! \class OPB< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  convection-diffusion operator
        (Laplacian u) + rho*(du / dx)                     
  on the unit squre [0,1]x[0,1] with zero Dirichlet boundary condition.                                                
*/
template <class ScalarType>
class OPB : public Anasazi::Operator<ScalarType>
{
private:
  ScalarType _rho;
  void tv(int nx, const ScalarType *x, ScalarType *y) const {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    ScalarType h = ONE / (ScalarType)(nx+1),
               h2 = h*h,
               dd = ((ScalarType)4.0) / h2,
               dl = -ONE / h2 - ((ScalarType)0.5)*_rho / h,
               du = -ONE / h2 + ((ScalarType)0.5)*_rho / h;
    
    
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
  
  OPB(const ScalarType rho) : _rho(rho) {}
  ~OPB() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    int n, nx;
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // deduce the size of the operator from the vector length...
    n = X.GetVecLength();
    // ... and the number of interior points in the discretization from that
    nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (nx*nx != n) return(Anasazi::Failed);
    
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



/******************************************************************************/
/*! \class OPC< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  convection-diffusion operator
          (d^2u/dx^2) + rho*(du/dx)                          
  on the interval [0,1] with zero Dirichlet boundary condition.  
*/
template <class ScalarType>
class OPC : public Anasazi::Operator<ScalarType>
{
private:
  ScalarType _rho;
public:
  
  OPC(const ScalarType rho) : _rho(rho) {}
  ~OPC() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    int n = X.GetVecLength();
    
    // Perform  Y <--- A*X 
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType TWO = (ScalarType)(2.0);
    ScalarType h  = ONE / (ScalarType)(n+1),
               s  = _rho*h / TWO,
               dd = TWO,
               dl = -ONE-s,
               du = -ONE+s;
    
    int j, p;
    for (p=0; p<nvecs; p++) {
      ScalarType *y = (*MyY)[p];
      const ScalarType *x = (*MyX)[p];
      j = 0;
      y[j] = dd*x[j] + du*x[j+1];
      for (j=1; j<n-1; j++) {
        y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
      }
      j = n-1;
      y[j] = dl*x[j-1] + dd*x[j];
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPD< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  convection-diffusion operator
          (d^2u/dx^2) + rho*(du/dx)                          
  on the interval [0,1] with zero Dirichlet boundary condition.  
  
  The operator applied is:
    OPD = inv[A-sigma*I]
*/
template <class ScalarType>
class OPD : public Anasazi::Operator<ScalarType>
{
private:
  int _n,_nx;
  ScalarType _rho, _sigma;
  std::vector<ScalarType> _dl, _dd, _du, _du2;
  std::vector<int> _ipiv;
  int _ferror;
  
public:
  
  OPD( const int n, const ScalarType rho, const ScalarType sigma) : _n(n), _rho(rho), _sigma(sigma) {
    
    typedef ScalarTraits<ScalarType> SCT;
    const ScalarType ONE = SCT::one();
    const ScalarType TWO = (ScalarType)(2.0)*ONE;
    LAPACK<int,ScalarType> lapack;
    
    _nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (_nx*_nx != n) {
      cout << "Argument 1 to OPD() was not a square number." << endl;
      _n = 100;
      _nx = 10;
    }
    
    /*----------------------------------------------------*
    | Construct C = A - SIGMA*I and factor C using LAPACK |
    | subroutine gttrf. The matrix A is chosen to be      |
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
    _dl.resize(_n-1,s1);
    _dd.resize(_n  ,s2);
    _du.resize(_n-1,s3);
    _du2.resize(_n-2);
    _ipiv.resize(_n);
  
    int _ferror;
    lapack.GTTRF(_n,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&_ferror);
    if (_ferror != 0) {
      cout << "Error in GTTRF in OPD()" << endl;
    }
  }
  ~OPD() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    LAPACK<int,ScalarType> lapack;
    
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
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    if (X.GetVecLength() != _n) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X = inv[A-SIGMA*I]*X using GTTRS
    int p;
    // set Y = X, as GTTRS operates in situ
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
    // call GTTRS multiple times (it takes multiple RHS, but MyMultiVec doesn't
    // use block storage)
    int ierr;
    for (p=0; p<nvecs; p++) {
      lapack.GTTRS('N',_n,1,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],(*MyY)[p],_n,&ierr);
      if (ierr != 0) return Anasazi::Failed;
    }
    
    return(Anasazi::Ok);
  }
};






/******************************************************************************/
/*! \class OPM< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of 2-D Laplacian.
*/
template <class ScalarType>
class OPM : public Anasazi::Operator<ScalarType>
{
private:
  void tv(int nx, const ScalarType *x, ScalarType *y) const {
    ScalarType dd = 4.0,
               dl = -ScalarTraits<ScalarType>::one(),
               du = -ScalarTraits<ScalarType>::one();
    
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
  
  OPM()  {}
  ~OPM() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    int n, nx;
  
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // deduce the size of the operator from the vector length...
    n = X.GetVecLength();
    // ... and the number of interior points in the discretization from that
    nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (nx*nx != n) return(Anasazi::Failed);
    
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
    
    // scale the vector by (1/h^2), where h is the mesh size
    ScalarType h2 = ONE/(ScalarType)((nx+1)*(nx+1));
    for (p=0; p<nvecs; p++) {
      blas.SCAL(n,ONE/h2,&(*MyY)[p][0],1);
    }
    
    return(Anasazi::Ok);
  }
  
};



/******************************************************************************/
/*! \class OPN< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  1-D discrete Laplacian on the interval [0,1] with zero Dirichlet b.c.
*/
template <class ScalarType>
class OPN : public Anasazi::Operator<ScalarType>
{
public:
  
  OPN() {}
  ~OPN() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType TWO = ((ScalarType)2.0)*ONE;
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int n = X.GetVecLength();
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X, a tridiagonal matrix multiply
    ScalarType dl = -ONE,
               dd = TWO,
               du = -ONE,
               h2 = ONE / ((ScalarType)(n+1)*(n+1));
    int p, j;
    for (p=0; p<nvecs; p++) {
      ScalarType *y = (*MyY)[p];
      const ScalarType *x = (*MyX)[p];
      j = 0;
      y[j] = dd*x[j] + du*x[j+1];
      for (j=1; j<n-1; j++) {
        y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
      }
      j = n-1;
      y[j] = dl*x[j-1] + dd*x[j];
      blas.SCAL(n, ONE/h2, y, 1);
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPO< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the central difference discretization of a 
  1-D Laplacian on the interval [0,1] with zero Dirichlet b.c.
  
  The operator applied is:
    OPO = inv[A-sigma*I]
  where A is as in OPN
*/
template <class ScalarType>
class OPO : public Anasazi::Operator<ScalarType>
{
private:
  int _n,_nx;
  ScalarType _sigma;
  std::vector<ScalarType> _dl, _dd, _du, _du2;
  std::vector<int> _ipiv;
  int _ferror;
  
public:
  
  OPO( const int n, const ScalarType sigma) : _n(n), _sigma(sigma) {
    
    typedef ScalarTraits<ScalarType> SCT;
    const ScalarType ONE = SCT::one();
    const ScalarType TWO = (ScalarType)(2.0)*ONE;
    LAPACK<int,ScalarType> lapack;
    
    _nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (_nx*_nx != n) {
      cout << "Argument 1 to OPO() was not a square number." << endl;
      _n = 100;
      _nx = 10;
    }
    
    /*----------------------------------------------------*
    | Construct C = A - SIGMA*I and factor C using LAPACK |
    | subroutine gttrf. The matrix A is chosen to be      |
    | the tridiagonal matrix derived from standard        |
    | central difference of the 1-d Laplacian             |
    \----------------------------------------------------*/
    ScalarType h2;
    h2 = ONE / (ScalarType)( (_n+1)*(_n+1) );
    _dd.resize(_n  , TWO/h2 - sigma );
    _dl.resize(_n-1, -ONE/h2 );
    _du.resize(_n-1, -ONE/h2 );
    _du2.resize(_n-2);
    _ipiv.resize(_n);
  
    int _ferror;
    lapack.GTTRF(_n,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&_ferror);
    if (_ferror != 0) {
      cout << "Error in GTTRF in OPO()" << endl;
    }
  }
  ~OPO() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    LAPACK<int,ScalarType> lapack;
    
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
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    if (X.GetVecLength() != _n) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X = inv[A-SIGMA*I]*X using GTTRS
    int p;
    // set Y = X, as GTTRS operates in situ
    MyY->MvAddMv( ONE, *MyX, ZERO, *MyX );
    // call GTTRS multiple times (it takes multiple RHS, but MyMultiVec doesn't
    // use block storage)
    int ierr;
    for (p=0; p<nvecs; p++) {
      lapack.GTTRS('N',_n,1,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],(*MyY)[p],_n,&ierr);
      if (ierr != 0) return Anasazi::Failed;
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPP< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the finite element discretization of a 
  1-D discrete Laplacian on the interval [0,1] with zero Dirichlet b.c.
  This is the A (stiffness) matrix.
*/
template <class ScalarType>
class OPP : public Anasazi::Operator<ScalarType>
{
public:
  
  OPP() {}
  ~OPP() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType TWO = ((ScalarType)2.0)*ONE;
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int n = X.GetVecLength();
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X, a tridiagonal matrix multiply
    ScalarType dl = -ONE,
               dd = TWO,
               du = -ONE,
               h  = ONE / ((ScalarType)(n+1));
    int p, j;
    for (p=0; p<nvecs; p++) {
      ScalarType *y = (*MyY)[p];
      const ScalarType *x = (*MyX)[p];
      j = 0;
      y[j] = dd*x[j] + du*x[j+1];
      for (j=1; j<n-1; j++) {
        y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
      }
      j = n-1;
      y[j] = dl*x[j-1] + dd*x[j];
      blas.SCAL(n, ONE/h, y, 1);
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPQ< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the finite element discretization of a 
  1-D discrete Laplacian on the interval [0,1] with zero Dirichlet b.c.
  This is the M (mass) matrix.
*/
template <class ScalarType>
class OPQ : public Anasazi::Operator<ScalarType>
{
public:
  
  OPQ() {}
  ~OPQ() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType FOUR = ((ScalarType)4.0)*ONE;
    const ScalarType SIX = ((ScalarType)6.0)*ONE;
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    if (MyX == 0) return Anasazi::Failed;
      
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    if (MyY == 0) return Anasazi::Failed;
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    
    int n = X.GetVecLength();
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X, a tridiagonal matrix multiply
    ScalarType dl = ONE,
               dd = FOUR,
               du = ONE,
               h  = ONE / ((ScalarType)(n+1)*SIX);
    int p, j;
    for (p=0; p<nvecs; p++) {
      ScalarType *y = (*MyY)[p];
      const ScalarType *x = (*MyX)[p];
      j = 0;
      y[j] = dd*x[j] + du*x[j+1];
      for (j=1; j<n-1; j++) {
        y[j] = dl*x[j-1] + dd*x[j] + du*x[j+1];
      }
      j = n-1;
      y[j] = dl*x[j-1] + dd*x[j];
      blas.SCAL(n, ONE/h, y, 1);
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPR< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the finite element discretization of a 
  1-D Laplacian on the interval [0,1] with zero Dirichlet b.c.
  
  The operator applied is:
    OPR = inv[M]*A
  where A is as in OPP and M is as in OPQ
*/
template <class ScalarType>
class OPR : public Anasazi::Operator<ScalarType>
{
private:
  int _n,_nx;
  std::vector<ScalarType> _dl, _dd, _du, _du2;
  std::vector<int> _ipiv;
  int _ferror;
  RefCountPtr< Anasazi::Operator<ScalarType> > _A;
  
public:
  
  OPR( const int n ) : _n(n) {
    
    typedef ScalarTraits<ScalarType> SCT;
    const ScalarType ONE = SCT::one();
    const ScalarType FOUR = ((ScalarType)4.0)*ONE;
    const ScalarType SIX = ((ScalarType)6.0)*ONE;
    LAPACK<int,ScalarType> lapack;
    
    // instantiate an A matrix 
    _A = rcp( new OPP<ScalarType>() );
    
    _nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (_nx*_nx != n) {
      cout << "Argument 1 to OPR() was not a square number." << endl;
      _n = 100;
      _nx = 10;
    }
    
    /*----------------------------------------------------*
    | Construct M and factor using LAPACK subroutine      |
    | gttrf. The matrix M is the tridiagonal matrix       |
    | derived from the finite element discretization      |
    | of the 1-d Laplacian.                               |
    \----------------------------------------------------*/
    ScalarType h, r1, r2;
    h = ONE / ((ScalarType)(_n+1));
    r1 = (FOUR / SIX) * h;
    r2 = (ONE / SIX) * h;
    _dd.resize(_n  , r1);
    _dl.resize(_n-1, r2);
    _du.resize(_n-1, r2);
    _du2.resize(_n-2);
    _ipiv.resize(_n);
  
    int _ferror;
    lapack.GTTRF(_n,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&_ferror);
    if (_ferror != 0) {
      cout << "Error in GTTRF in OPR()" << endl;
    }
  }
  ~OPR() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    LAPACK<int,ScalarType> lapack;
    
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
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    if (X.GetVecLength() != _n) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X = inv[A-SIGMA*I]*X using GTTRS
    int p;
    // set Y = A*X, as GTTRS operates in situ
    if ( _A->Apply(*MyX,*MyY) != Anasazi::Ok ) return Anasazi::Failed;
    // now, perform inv[M]*Y = inv[M]*X
    // call GTTRS multiple times (it takes multiple RHS, but MyMultiVec doesn't
    // use block storage)
    int ierr;
    for (p=0; p<nvecs; p++) {
      lapack.GTTRS('N',_n,1,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],(*MyY)[p],_n,&ierr);
      if (ierr != 0) return Anasazi::Failed;
    }
    
    return(Anasazi::Ok);
  }
};



/******************************************************************************/
/*! \class OPS< ScalarType >
  \brief Implementation of Anasazi::MultiVector< ScalarType > for the
  application of the finite element discretization of a 1-D Laplacian
  on the interval [0,1] with zero Dirichlet b.c.
  
  The operator applied is:
    OPS = inv[A-sigma*M]*M
  where A is as in OPP and M is as in OPQ
*/
template <class ScalarType>
class OPS : public Anasazi::Operator<ScalarType>
{
private:
  int _n,_nx;
  ScalarType _sigma;
  std::vector<ScalarType> _dl, _dd, _du, _du2;
  std::vector<int> _ipiv;
  int _ferror;
  RefCountPtr< Anasazi::Operator<ScalarType> > _M;
  
public:
  
  OPS( const int n, const ScalarType sigma) : _n(n), _sigma(sigma) {
    
    typedef ScalarTraits<ScalarType> SCT;
    const ScalarType ONE = SCT::one();
    const ScalarType TWO = (ScalarType)(2.0)*ONE;
    const ScalarType FOUR = (ScalarType)(4.0)*ONE;
    const ScalarType SIX = (ScalarType)(6.0)*ONE;
    LAPACK<int,ScalarType> lapack;

    _M = rcp( new OPQ<ScalarType>() );
    
    _nx = ScalarTraits<int>::squareroot(n);
    // return an error if the vector length isn't a square number
    if (_nx*_nx != n) {
      cout << "Argument 1 to OPS() was not a square number." << endl;
      _n = 100;
      _nx = 10;
    }
    
    /*----------------------------------------------------*
    | Construct C = A - SIGMA*I and factor C using LAPACK |
    | subroutine gttrf.                                   |
    \----------------------------------------------------*/
    ScalarType h,r1,r2;
    h = ONE / (ScalarType)(_n+1);
    r1 = (FOUR / SIX) * h;
    r2 = (ONE / SIX) * h;
    
    _dl.resize(_n-1, -ONE/h - sigma*r2 );
    _dd.resize(_n  ,  TWO/h - sigma*r1 );
    _du.resize(_n-1, -ONE/h - sigma*r2 );
    _du2.resize(_n-2);
    _ipiv.resize(_n);
  
    int _ferror;
    lapack.GTTRF(_n,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],&_ferror);
    if (_ferror != 0) {
      cout << "Error in GTTRF in OPS()" << endl;
    }
    
    /*
    typename std::vector<ScalarType>::iterator it;
    typename std::vector<int>::iterator it2;
    cout << "_dl" << endl;
    for (it=_dl.begin(); it != _dl.end(); it++) cout << *it << endl;
    cout << "_dd" << endl;
    for (it=_dd.begin(); it != _dd.end(); it++) cout << *it << endl;
    cout << "_du" << endl;
    for (it=_du.begin(); it != _du.end(); it++) cout << *it << endl;
    cout << "_du2" << endl;
    for (it=_du2.begin(); it != _du2.end(); it++) cout << *it << endl;
    cout << "ipiv" << endl;
    for (it2=_ipiv.begin(); it2 != _ipiv.end(); it2++) cout << *it2 << endl;
    */
    
  }
  ~OPS() {}
  
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const ScalarType ONE = ScalarTraits<ScalarType>::one();
    const ScalarType ZERO = ScalarTraits<ScalarType>::zero();
    BLAS<int,ScalarType> blas;
    LAPACK<int,ScalarType> lapack;
    
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
      
    if (X.GetNumberVecs() != Y.GetNumberVecs()) return Anasazi::Failed;
    if (X.GetVecLength() != Y.GetVecLength()) return Anasazi::Failed;
    if (X.GetVecLength() != _n) return Anasazi::Failed;
    
    int nvecs = X.GetNumberVecs();
    
    // Perform  Y <--- OP*X = inv[A-SIGMA*I]*X using GTTRS
    int p;
    // set Y = M*X, as GTTRS operates in situ
    _M->Apply( *MyX, *MyY );
    // set Y = inv[A-sigma*M]*Y = inv[A-sigma*M]*M*X
    // call GTTRS multiple times (it takes multiple RHS, but MyMultiVec doesn't
    // use block storage)
    int ierr;
    for (p=0; p<nvecs; p++) {
      lapack.GTTRS('N',_n,1,&_dl[0],&_dd[0],&_du[0],&_du2[0],&_ipiv[0],(*MyY)[p],_n,&ierr);
      if (ierr != 0) return Anasazi::Failed;
    }
    
    return(Anasazi::Ok);
  }
};


template <class ScalarType>
class ARPACK_Example {
  public:
    virtual void xformeval(std::vector<ScalarType> &) const = 0;
    virtual RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const = 0;
    virtual RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const = 0;
    virtual RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const = 0;
    virtual bool        isHerm() const = 0;
    virtual std::string getSort() const = 0;
};


template <class ScalarType>
class ARPACK_NDRV1 : public ARPACK_Example<ScalarType> {
  private:
    ScalarType _rho;
  public:
    ARPACK_NDRV1(ScalarType rho = ScalarTraits<ScalarType>::zero()) : _rho(rho) {}
    
    void xformeval(std::vector<ScalarType> &vals) const {}
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPB<ScalarType>(_rho)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPB<ScalarType>(_rho)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPA<ScalarType>()); }
    bool isHerm() const {return false;}
    std::string getSort() const {return string("SM");}
};

template <class ScalarType>
class ARPACK_NDRV2 : public ARPACK_Example<ScalarType> {
  private:
    ScalarType _rho, _sigma;
    int _n;
  public:
    ARPACK_NDRV2(int n,
                 ScalarType rho = (ScalarType)(1.0e+1),
                 ScalarType sigma = (ScalarType)(1.0) ) : _n(n), _rho(rho),_sigma(sigma) {}
    
    void xformeval(std::vector<ScalarType> &vals) const { 
      typename std::vector<ScalarType>::iterator i;
      const ScalarType ONE = ScalarTraits<ScalarType>::one();
      for (i=vals.begin(); i!=vals.end(); i++) {
        *i = ONE / *i + _sigma;
      }
    }
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPD<ScalarType>(_n,_rho,_sigma)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPC<ScalarType>(_rho)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPA<ScalarType>()); }
    bool isHerm() const {return false;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_NDRV3 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return false;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_NDRV4 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return false;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_NDRV5 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return false;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_NDRV6 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return false;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_SDRV1 : public ARPACK_Example<ScalarType> {
  private:
  public:
    void xformeval(std::vector<ScalarType> &vals) const {}
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPM<ScalarType>()); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPM<ScalarType>()); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPA<ScalarType>()); }
    bool isHerm() const {return true;}
    std::string getSort() const {return string("SM");}
};

template <class ScalarType>
class ARPACK_SDRV2 : public ARPACK_Example<ScalarType> {
  private:
    ScalarType _sigma;
    int _n;
  public:
    ARPACK_SDRV2(int n, ScalarType sigma = ScalarTraits<ScalarType>::zero()) 
        : _n(n), _sigma(sigma) {}
    void xformeval(std::vector<ScalarType> &vals) const {
      typename std::vector<ScalarType>::iterator i;
      const ScalarType ONE = ScalarTraits<ScalarType>::one();
      for (i=vals.begin(); i!=vals.end(); i++) {
        *i = ONE / *i + _sigma;
      }
    }
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPO<ScalarType>(_n,_sigma)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPN<ScalarType>()); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPA<ScalarType>()); }
    bool isHerm() const {return true;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_SDRV3 : public ARPACK_Example<ScalarType> {
  private:
    int _n;
  public:
    ARPACK_SDRV3(int n) : _n(n) {}
    void xformeval(std::vector<ScalarType> &vals) const {}
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPR<ScalarType>(_n)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPP<ScalarType>()); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPQ<ScalarType>()); }
    bool isHerm() const {return true;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_SDRV4 : public ARPACK_Example<ScalarType> {
  private:
    int _n;
    ScalarType _sigma;
  public:
    ARPACK_SDRV4(int n, ScalarType sigma = ScalarTraits<ScalarType>::zero()) 
        : _n(n), _sigma(sigma) {}
    void xformeval(std::vector<ScalarType> &vals) const {
      typename std::vector<ScalarType>::iterator i;
      const ScalarType ONE = ScalarTraits<ScalarType>::one();
      for (i=vals.begin(); i!=vals.end(); i++) {
        *i = ONE / *i + _sigma;
      }
    }
    RefCountPtr< Anasazi::Operator<ScalarType> > getOp() const { return rcp(new OPS<ScalarType>(_n,_sigma)); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getA()  const { return rcp(new OPP<ScalarType>()); }
    RefCountPtr< Anasazi::Operator<ScalarType> > getM()  const { return rcp(new OPQ<ScalarType>()); }
    bool isHerm() const {return true;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_SDRV5 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return true;}
    std::string getSort() const {return string("LM");}
};

template <class ScalarType>
class ARPACK_SDRV6 : public ARPACK_Example<ScalarType> {
  private:
  public:
    bool isHerm() const {return true;}
    std::string getSort() const {return string("LM");}
};


/*

ARPACK_NDRV3
OPG
OPE
OPF

ARPACK_NDRV4
OPH
OPE
OPF

ARPACK_NDRV5
OPK
OPI
OPJ

ARPACK_NDRV6
OPL
OPI
OPJ

ARPACK_SDRV1
OPM
OPM
OPA

ARPACK_SDRV3
OPR
OPP
OPQ

ARPACK_SDRV4
OPS
OPP
OPQ

ARPACK_SDRV5
OPT
OPP
OPP

ARPACK_SDRV6
OPU
OPP
OPQ
*/

template <class ScalarType>
RefCountPtr< ARPACK_Example<ScalarType> > GetARPACKExample(const std::string &drivername, int dim) {
  RefCountPtr< ARPACK_Example<ScalarType> > nullptr;

  std::string dncopy(drivername);
  // if they sent the full driver name, remove the scalar type [sdcz]
  if (dncopy.length() == 6 && dncopy.find_first_of("SDCZsdcz") == 0) {
    dncopy = dncopy.substr(1);
  }

  if (dncopy == "ndrv1" || dncopy == "NDRV1") {
    return rcp( new ARPACK_NDRV1<ScalarType>() );
  }
  else if (dncopy == "ndrv2" || dncopy == "NDRV2") {
    return rcp( new ARPACK_NDRV2<ScalarType>(dim) );
  }
  else if (dncopy == "sdrv1" || dncopy == "SDRV1") {
    return rcp( new ARPACK_SDRV1<ScalarType>() );
  }
  else if (dncopy == "sdrv2" || dncopy == "SDRV2") {
    return rcp( new ARPACK_SDRV2<ScalarType>(dim) );
  }
  else if (dncopy == "sdrv3" || dncopy == "SDRV3") {
    return rcp( new ARPACK_SDRV3<ScalarType>(dim) );
  }
  else if (dncopy == "sdrv4" || dncopy == "SDRV4") {
    return rcp( new ARPACK_SDRV4<ScalarType>(dim) );
  }
  return nullptr;
}


#endif //ARPACK_OPERATORS_HPP
