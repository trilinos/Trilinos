
#ifndef TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP
#define TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

#include <Tpetra_CrsMatrix.hpp>

#define TIFPACK_INST(CLASSNAME,S,LO,GO) \
  template class CLASSNAME<Tpetra::CrsMatrix<S,LO,GO, \
                 Kokkos::DefaultNode::DefaultNodeType, \
                 Kokkos::DefaultKernels<S,LO,Kokkos::DefaultNode::DefaultNodeType>::SparseOps> >
  
#define TIFPACK_CLASS_CrsMatrix_float_int_int_defaultNode_defaultOps(CLASSNAME) \
  TIFPACK_INST(CLASSNAME,float,int,int)

#define TIFPACK_CLASS_CrsMatrix_float_short_int_defaultNode_defaultOps(CLASSNAME) \
  TIFPACK_INST(CLASSNAME,float,short,int)

#define TIFPACK_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME) \
  TIFPACK_INST(CLASSNAME,double,int,int)

#define TIFPACK_INSTANT_CRSMATRIX_FLOAT_DOUBLE_DEFAULTS(CLASSNAME) \
  TIFPACK_CLASS_CrsMatrix_double_int_int_defaultNode_defaultOps(CLASSNAME)

#define TIFPACK_INSTANT_CRSMATRIX_COMPLEX_DEFAULTS(CLASSNAME) \
  TIFPACK_INST(CLASSNAME,std::complex<double>,int,int) \
  TIFPACK_INST(CLASSNAME,std::complex<float>,int,int) 

#endif

