#ifndef ARPACK_OPS_HPP
#define ARPACK_OPS_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziOperator.hpp"
#include "MyMultiVec.hpp"

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
 * \author Chris Baker (SNL/1414,FSU/CSIT) and Heidi Thornquist (SNL/1437)
 *
 * \date Last modified on 12-Dec-05
 */


template <class ScalarType>
struct ARPACKEx1
{
  typedef   OPE<ScalarType>   OP;
  typedef   OPA<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx2
{
  typedef   OPF<ScalarType>   OP;
  typedef   OPA<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx3
{
  typedef   OPN<ScalarType>   OP;
  typedef   OPM<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx4
{
  typedef   OPO<ScalarType>   OP;
  typedef   OPM<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx5
{
  typedef   OPC<ScalarType>   OP;
  typedef   OPB<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx6
{
  typedef   OPD<ScalarType>   OP;
  typedef   OPB<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx7
{
  typedef   OPG<ScalarType>   OP;
  typedef   OPA<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx8
{
  typedef   OPI<ScalarType>   OP;
  typedef   OPA<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx9
{
  typedef   OPR<ScalarType>   OP;
  typedef   OPQ<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx10
{
  typedef   OPS<ScalarType>   OP;
  typedef   OPQ<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx11
{
  typedef   OPU<ScalarType>   OP;
  typedef   OPP<ScalarType>   B;
}

template <class ScalarType>
struct ARPACKEx12
{
  typedef   OPT<ScalarType>   OP;
  typedef   OPQ<ScalarType>   B;
}



template <class ScalarType>
class OPA : public Anasazi::Operator<ScalarType>
{
  
public:
  
  OPA(const int numRows) :
    _numRows(numRows)
  {}
  
  //! Dtor
  ~OPA()
  {}
  
  //! Applies the operator
  Anasazi::ReturnType Apply(const Anasazi::MultiVec<ScalarType>& X, 
                                  Anasazi::MultiVec<ScalarType>& Y ) const
  {
    const MyMultiVec<ScalarType>* MyX;
    MyX = dynamic_cast<const MyMultiVec<ScalarType>*>(&X); 
    assert (MyX != 0);
    
    MyMultiVec<ScalarType>* MyY;
    MyY = dynamic_cast<MyMultiVec<ScalarType>*>(&Y); 
    assert (MyY != 0);
    
    assert (X.GetNumberVecs() == Y.GetNumberVecs());
    assert (X.GetVecLength() == Y.GetVecLength());
   
    return(Anasazi::Ok);
  }
  
private:
  //! Number of rows and columns
  int NumRows_;
};

#endif //ARPACK_OPS_HPP
