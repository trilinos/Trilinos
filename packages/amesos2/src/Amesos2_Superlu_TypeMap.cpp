/**
 * \file   Amesos2_Superlu_TypeMap.cpp
 * \author Eric Bavier <etbavie@sandia.gov>
 * \date   Fri Jul 22 11:25:30 2011
 * 
 * \brief  Definitions for SuperLU TypeMap.
 */

#include "Amesos2_Superlu_TypeMap.hpp"

namespace Amesos2 {
  
  SLU::Dtype_t TypeMap<Superlu,float>::dtype = SLU::SLU_S;

  SLU::Dtype_t TypeMap<Superlu,double>::dtype = SLU::SLU_D;

#ifdef HAVE_TEUCHOS_COMPLEX
  SLU::Dtype_t TypeMap<Superlu,std::complex<float> >::dtype = SLU::SLU_C;

  SLU::Dtype_t TypeMap<Superlu,std::complex<double> >::dtype = SLU::SLU_Z;

  SLU::Dtype_t TypeMap<Superlu,SLU::C::complex>::dtype = SLU::SLU_C;

  SLU::Dtype_t TypeMap<Superlu,SLU::Z::doublecomplex>::dtype = SLU::SLU_Z;
#endif
  
}

#ifdef HAVE_TEUCHOS_COMPLEX
namespace std {
  ostream& operator<<(ostream& out, const SLU::Z::doublecomplex z){
    return (out << "(" << z.r << "," << z.i << ")");
  }

  ostream& operator<<(ostream& out, const SLU::C::complex c){
    return (out << "(" << c.r << "," << c.i << ")");
  }
}
#endif
