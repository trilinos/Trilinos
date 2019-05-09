#include "Intrepid_ArrayTools.hpp"
#include <iostream>
#ifdef INTREPID_OLD_KOKKOS_CODE
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <impl/Kokkos_Timer.hpp>
#endif
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_ScalarTraits.hpp"
#ifdef KOKKOS_ENABLE_CUDA
#include "cublas_v2.h"
#endif
#include <Teuchos_BLAS.hpp>
#include <Teuchos_as.hpp>

namespace KokkosDenseMat{
	template<typename Scalar, typename DeviceType,typename Layout,int rank>
	struct MultiGemm;
	
template<typename Scalar>
	struct MultiGemm<Scalar,void,void,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, int m, int n, int k, Scalar alpha,
          Scalar* A, Scalar* B,
          Scalar beta, Scalar* C){
	Teuchos::BLAS<int,Scalar>blas;

	blas.GEMM(transA, transB, n, m, k, alpha,
                   B, n,
                   A, k,
                   beta, C, n);
					   
					   
	}
	
};

template<typename Scalar>
	struct MultiGemm<Scalar,void,void,1>{
	static Scalar DDOT(int n, Scalar* A, int incA, Scalar* B, int incB){
	Teuchos::BLAS<int,Scalar>blas;
	Scalar result = blas.DOT(n, A, incA, B, incB);
	return result;				   
	}
	
};
#ifdef INTREPID_OLD_KOKKOS_CODE	

template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
          Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> A, Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> C){
	Teuchos::BLAS<int,Scalar>blas;
	    const int m = static_cast<int> (C.dimension_0 ()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());

	blas.GEMM(transA, transB, m, n, k, alpha,
                   A.ptr_on_device(), k,
                   B.ptr_on_device(), n,
                   beta, C.ptr_on_device(), n);
					   
					   
	}
	
};

template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,2>{
	static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
          Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> A, Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> C){
	Teuchos::BLAS<int,Scalar>blas;
	    const int m = static_cast<int> (C.dimension_0 ()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
   
	blas.GEMM(transB, transA, n, m, k, alpha,
                   B.ptr_on_device(), n,
                   A.ptr_on_device(), k,
                   beta, C.ptr_on_device(), n);
					   
					   
	}
	
};



template<class Scalar>
struct blasOpenMPBatchLeft {
typedef Kokkos::View<Scalar***, Kokkos::LayoutLeft> left_type;
typedef Kokkos::View<Scalar***, Kokkos::LayoutRight> right_type;
typedef Kokkos::View<Scalar**, Kokkos::LayoutStride> left_subtype;
 left_type A;
 left_type B;
 left_type C;
  
  

  int msize,nsize,ksize;
  Teuchos::ETransp transA;
  Teuchos::ETransp transB;
  Scalar alpha;
  Scalar beta;
  blasOpenMPBatchLeft (left_type a_, left_type b_, left_type c_, const int m_, const int n_, const int k_,
  Teuchos::ETransp transA_, Teuchos::ETransp transB_, Scalar alpha_, Scalar beta_) :
    A(a_), B(b_), C(c_), msize(m_), nsize(n_), ksize(k_), transA(transA_), transB(transB_),
    alpha(alpha_), beta(beta_)
  {}
 
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {

	 left_subtype subA=subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subB=subview(B, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subC=subview(C, i, Kokkos::ALL(), Kokkos::ALL());

//#pragma loop(ivdep)

for (size_t i = 0; i < msize; i++) {
  
    for (size_t k = 0; k < ksize; k++) {

        const Scalar r = subA(i,k);

        for (size_t j = 0; j < nsize; j++) {

		subC(i,j)=beta*subC(i,j)+alpha*r*subB(k,j);

        }

    }

} 
                        
  }
};

template<class Scalar>
struct blasOpenMPBatchRight {
typedef Kokkos::View<Scalar***, Kokkos::LayoutLeft> left_type;
typedef Kokkos::View<Scalar***, Kokkos::LayoutRight> right_type;
typedef Kokkos::View<Scalar**, Kokkos::LayoutStride> left_subtype;

 Teuchos::BLAS<int,Scalar>blas;
 right_type A;
 right_type B;
 right_type C;
  
  

  int msize,nsize,ksize;
  Teuchos::ETransp transA;
  Teuchos::ETransp transB;
  Scalar alpha;
  Scalar beta;
  blasOpenMPBatchRight (right_type a_, right_type b_, right_type c_, const int m_, const int n_, const int k_,
  Teuchos::ETransp transA_, Teuchos::ETransp transB_, Scalar alpha_, Scalar beta_) :
    A(a_), B(b_), C(c_), msize(m_), nsize(n_), ksize(k_), transA(transA_), transB(transB_),
    alpha(alpha_), beta(beta_)
  {}
 
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {

	 left_subtype subA=subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subB=subview(B, i, Kokkos::ALL(), Kokkos::ALL());
	 left_subtype subC=subview(C, i, Kokkos::ALL(), Kokkos::ALL());

//#pragma loop(ivdep)

for (std::size_t i = 0; i != static_cast<std::size_t>(msize); ++i) {
  
    for (std::size_t k = 0; k != static_cast<std::size_t>(ksize); ++k) {

        const Scalar r = subA(i,k);

        for (std::size_t j = 0; j != static_cast<std::size_t>(nsize); ++j) {

		subC(i,j)=beta*subC(i,j)+alpha*r*subB(k,j);

        }

    }

} 
                        
  }
};


template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
           Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> A,  Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar***,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
   //     printf("m:%d,n:%d,k:%d",m,n,k);
	Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchLeft<Scalar>(A,B,C,m,n,k,transA,transB,alpha,beta));
	}
	
};
template<typename Scalar>
	struct MultiGemm<Scalar,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, Scalar alpha,
           Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> A,  Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> B,
          Scalar beta, Kokkos::View<Scalar***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
Teuchos::BLAS<int,Scalar>blas;
Kokkos::parallel_for(C.dimension_0(),KOKKOS_LAMBDA (const size_t i) {
        blas.GEMM(transB, transA, n, m, k, alpha,
                   &B(i,0,0), n,
                   &A(i,0,0), k,
                   beta, &C(i,0,0), n);
});
	}
	
};
#endif

// The Following stuff doesn't compile at all with Cuda. And it can't for example Scalar is not defined because its explicit specialisations.
#ifdef KOKKOS_ENABLE_CUDA
template <>	
	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutLeft,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
stat = cublasCreate(&handle);
stat = cublasDgemm(handle,transposeA, transposeB, m, n, k, &alpha, A.ptr_on_device(), k, B.ptr_on_device(), n, &beta, C.ptr_on_device(), n);
	}
	
};
template <>
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutLeft,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);
  stat=cublasSgemm(handle,transposeA, transposeB, m, n, k, &alpha, A.ptr_on_device(), k, B.ptr_on_device(), n, &beta, C.ptr_on_device(), n);
	}
	
};

template <>
	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutRight,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          double beta, Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle); 
  stat= cublasDgemm(handle,transposeB, transposeA, n, m, k, &alpha, B.ptr_on_device(), n, A.ptr_on_device(), k, &beta, C.ptr_on_device(), n);
	}
	
};

template <>
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutRight,2>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> B,
          float beta, Kokkos::View<float**,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_0()),
        n = static_cast<int> (C.dimension_1 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_1 () : A.dimension_0 ());
cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
 stat = cublasCreate(&handle);  
 stat= cublasSgemm(handle,transposeB, transposeA, n, m, k, &alpha, B.ptr_on_device(), n, A.ptr_on_device(), k, &beta, C.ptr_on_device(), n);

	}
	
};

template <>
	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          double beta, Kokkos::View<double***,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
      
 Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchLeft<double>(A,B,C,m,n,k,transA,transB,alpha,beta));
	}
	
};

template <>	
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutLeft,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> A,  Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> B,
          float beta, Kokkos::View<float***,Kokkos::LayoutLeft,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());

 Kokkos::parallel_for(C.dimension(0),blasOpenMPBatchLeft<float>(A,B,C,m,n,k,transA,transB,alpha,beta));


	}
	
};

template <>
	struct MultiGemm<double,Kokkos::Cuda,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, double alpha,
           Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> B,
          double beta, Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());

cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);
size_t batchCount=C.dimension_0();
double ** Aptrs = (double ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(double * ));
double ** Bptrs = (double ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(double * ));
double ** Cptrs = (double ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(double * ));
for(size_t i=0;i<batchCount;i++){
Aptrs[i]=&A(i,0,0);
Bptrs[i]=&B(i,0,0);
Cptrs[i]=&C(i,0,0);
}
stat=cublasDgemmBatched(handle,transposeB,  transposeA,
                                  n, m, k,
                                  &alpha,
                                  (const double**)Bptrs, n,
                                  (const double**)Aptrs, k,
                                  &beta,
                                  Cptrs, n, batchCount);	
Kokkos::kokkos_free(Aptrs); 
Kokkos::kokkos_free(Bptrs);
Kokkos::kokkos_free(Cptrs);
}
	
};



template <>
	struct MultiGemm<float,Kokkos::Cuda,Kokkos::LayoutRight,3>{
		static void GEMM(Teuchos::ETransp transA, Teuchos::ETransp transB, float alpha,
           Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> A,  Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> B,
          float beta, Kokkos::View<float***,Kokkos::LayoutRight,Kokkos::Cuda> C){
		const int m = static_cast<int> (C.dimension_1()),
        n = static_cast<int> (C.dimension_2 ()),
        k = (transA == Teuchos::NO_TRANS ? A.dimension_2 () : A.dimension_1 ());
cublasOperation_t transposeA;
cublasOperation_t transposeB;
if(transA==Teuchos::NO_TRANS){
transposeA=CUBLAS_OP_N;
}else{
transposeA=CUBLAS_OP_T;
}
if(transB==Teuchos::NO_TRANS){
transposeB=CUBLAS_OP_N;
}else{
transposeB=CUBLAS_OP_T;
}

  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);
size_t batchCount=C.dimension_0();
float ** Aptrs = (float ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(float * ));
float ** Bptrs = (float ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(float * ));
float ** Cptrs = (float ** ) Kokkos::kokkos_malloc <Kokkos::CudaUVMSpace >(batchCount * sizeof(float * ));
for(size_t i=0;i<batchCount;i++){
Aptrs[i]=&A(i,0,0);
Bptrs[i]=&B(i,0,0);
Cptrs[i]=&C(i,0,0);
}
stat=cublasSgemmBatched(handle,transposeB,  transposeA,
                                  n, m, k,
                                  &alpha,
                                  (const float**)Bptrs, n,
                                  (const float**)Aptrs, k,
                                  &beta,
                                  Cptrs, n, batchCount);

Kokkos::kokkos_free(Aptrs);
Kokkos::kokkos_free(Bptrs);
Kokkos::kokkos_free(Cptrs);
}
	
};
#endif

}





int main(){
#ifdef INTREPID_OLD_KOKKOS_CODE

   Kokkos::initialize();
  //initialize viewsto random values
   {
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> inputview1("X",5,5);
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> inputview2("Y",5,5);
	Kokkos::View<double**,Kokkos::LayoutLeft,Kokkos::DefaultExecutionSpace> outputview2("Z",5,5);

  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutLeft,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
    /* 
std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview1(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<outputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}*/
}// scope 1
{// scope 2 

	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview1("X",5,5);
	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview2("Y",5,5);
	Kokkos::View<double**,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> outputview2("Z",5,5);

  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);	
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
/*std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview1(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<inputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<outputview2(i,j)<<",";
		
	}
	std::cout <<std::endl;
}*/
}
//begin scope 4
{
Intrepid::FieldContainer<double>testcontainerA(5,5);
	for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		testcontainerA(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
Intrepid::FieldContainer<double>testcontainerB(5,5);
	for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		testcontainerB(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
	Intrepid::FieldContainer<double>testcontainerC(5,5);
	KokkosDenseMat::MultiGemm<double,void,void,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,5,5,5,1.0,&testcontainerA[0],&testcontainerB[0],1,&testcontainerC[0]);
	/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerA(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerB(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<5;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerC(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<std::endl;	
*/
}//end scope 4
{
	Intrepid::FieldContainer<double>testcontainerA(2,5,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		testcontainerA(i,j,k)=Teuchos::ScalarTraits<double>::random();
	     }
	}
}
Intrepid::FieldContainer<double>testcontainerB(2,5,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		testcontainerB(i,j,k)=Teuchos::ScalarTraits<double>::random();
		}
	}
}
	Intrepid::FieldContainer<double>testcontainerC(2,5,5);	
for(int i=0;i<2;i++){
	KokkosDenseMat::MultiGemm<double,void,void,2>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,5,5,5,1.0,&testcontainerA[i*5*5],&testcontainerB[i*5*5],1,&testcontainerC[i*5*5]);
}	
/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerA(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerB(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
	std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<testcontainerC(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}	*/
	
}//end scope5

//scope6
{
	Intrepid::FieldContainer<double>testcontainerA(5);
	for(int i=0;i<5;i++){
		testcontainerA(i)=Teuchos::ScalarTraits<double>::random();
}
Intrepid::FieldContainer<double>testcontainerB(5);
	for(int i=0;i<5;i++){
		testcontainerB(i)=Teuchos::ScalarTraits<double>::random();
}
/*double result=	KokkosDenseMat::MultiGemm<double,void,void,1>::DDOT(5, &testcontainerA[0], 1, &testcontainerB[0], 1);

std::cout <<"A: "<<std::endl;
	for(int i=0;i<5;i++){
		std::cout <<testcontainerA(i)<<",";
}
std::cout <<std::endl;

std::cout <<"B: "<<std::endl;
	for(int i=0;i<5;i++){
		std::cout <<testcontainerB(i)<<std::endl;
}
std::cout <<std::endl;

std::cout <<"result: "<<result<<std::endl;

*/
}//end scope6

//scope 7
{
Intrepid::FieldContainer<double>testcontainerA(2,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		testcontainerA(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}
Intrepid::FieldContainer<double>testcontainerB(2,5);
	for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		testcontainerB(i,j)=Teuchos::ScalarTraits<double>::random();
	}
}

Intrepid::FieldContainer<double>testcontainerC(2);
for(int i=0;i<2;i++){
	testcontainerC(i)=KokkosDenseMat::MultiGemm<double,void,void,1>::DDOT(5, &testcontainerA[i*5], 1, &testcontainerB[i*5], 1);
}
/*
std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerA(i,j)<<",";
		
	}
	std::cout <<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		std::cout <<testcontainerB(i,j)<<std::endl;
		
	}
	std::cout <<std::endl;
}
std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
		std::cout <<testcontainerC(i)<<",";		
	std::cout <<std::endl;
}*/
}	
{//end scope 8

    Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview1("X",1000,5,6);
	Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> inputview2("Y",1000,6,5);
	Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::DefaultExecutionSpace> outputview2("Z",1000,5,5);
  //fill with random numbers
	Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
    Kokkos::fill_random(inputview1,rand_pool64,100);
    Kokkos::fill_random(inputview2,rand_pool64,100);
    KokkosDenseMat::MultiGemm<double,Kokkos::DefaultExecutionSpace,Kokkos::LayoutRight,3>::GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,inputview1,inputview2,1,outputview2);
/*	std::cout <<"A:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<6;k++){
		std::cout <<inputview1(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
std::cout <<"B:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<6;j++){
		for (int k=0;k<5;k++){
		std::cout <<inputview2(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}
	std::cout <<"C:"<<std::endl;
for(int i=0;i<2;i++){
	for (int j=0;j<5;j++){
		for (int k=0;k<5;k++){
		std::cout <<outputview2(i,j,k)<<",";
	     }
	     std::cout <<std::endl;
	}
	std::cout <<std::endl<<std::endl;
}*/

}

   Kokkos::finalize();
  
#endif 
  
	
	return 0;
}
