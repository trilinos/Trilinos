
#ifndef TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP
#define TPETRA_EXPLICITINSTANTIATIONHELPERS_HPP

#include <Tpetra_CrsMatrix.hpp>

#define TIFPACK_CLASS_CrsMatrix_float_int_int_defaultNode_defaultSM_defaultSS(CLASSNAME) \
  template class CLASSNAME<Tpetra::CrsMatrix<float,int,int, \
                  Kokkos::DefaultNode::DefaultNodeType, \
                  Kokkos::DefaultKernels<float,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps > >;

#define TIFPACK_CLASS_CrsMatrix_float_short_int_defaultNode_defaultSM_defaultSS(CLASSNAME) \
  template class CLASSNAME<Tpetra::CrsMatrix<float,short,int, \
                  Kokkos::DefaultNode::DefaultNodeType, \
                  Kokkos::DefaultKernels<float,short,Kokkos::DefaultNode::DefaultNodeType>::SparseOps > >;

#define TIFPACK_CLASS_CrsMatrix_double_int_int_defaultNode_defaultSM_defaultSS(CLASSNAME) \
  template class CLASSNAME<Tpetra::CrsMatrix<double,int,int, \
                  Kokkos::DefaultNode::DefaultNodeType, \
                  Kokkos::DefaultKernels<double,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps > >;

#define TIFPACK_INSTANT_CRSMATRIX_FLOAT_DOUBLE_DEFAULTS(CLASSNAME) \
  TIFPACK_CLASS_CrsMatrix_double_int_int_defaultNode_defaultSM_defaultSS(CLASSNAME)

//  TIFPACK_CLASS_CrsMatrix_float_int_int_defaultNode_defaultSM_defaultSS(CLASSNAME)
//  TIFPACK_CLASS_CrsMatrix_float_short_int_defaultNode_defaultSM_defaultSS(CLASSNAME)

#endif

