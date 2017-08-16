/*
// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Kokkos_Core.hpp"
#include "Sacado.hpp"
#include <impl/Kokkos_Timer.hpp>

#include <random>
#include <time.h>
#include <stdlib.h>
#include "Intrepid2_FieldContainer.hpp"
#include <Kokkos_Random.hpp>
#if 0
namespace Intrepid{
class none{};
class FieldContainer_Kokkos_Ptr;
template <class Scalar,class ScalarPointer=void,class MemoryLayout=void,class ExecutionSpace=void>
class FieldContainer_Kokkos;
template <class Scalar>
class FieldContainer_Kokkos<Scalar,void,void,void>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};

Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();
Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num)const{return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0){delete[] containerMemory;}
}
template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar>
FieldContainer_Kokkos<Scalar>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}
#ifdef KOKKOS_HAVE_SERIAL
template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
unsigned int count_=1;
bool intrepidManaged=false;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();
FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>& InContainer);

Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();
FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>& InContainer);


Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}



template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
#endif

#ifdef KOKKOS_HAVE_OPENMP
template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();

FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>& InContainer);

Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}




template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}


template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();


FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>& InContainer);

Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}




template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::OpenMP>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
#endif

#ifdef KOKKOS_HAVE_PTHREAD

template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();


FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>& InContainer);

Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}



template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}


template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
Scalar* containerMemory;
index_type rankValue=0;
index_type sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();



FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>& InContainer);
Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){delete[] containerMemory;}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
}




template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();


}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

#endif

#ifdef KOKKOS_HAVE_CUDA

template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
bool intrepidManaged=false;
Scalar* containerMemory;
unsigned int count_=1;
index_type rankValue=0;
index_type sizeValue=0;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>& InContainer);

FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();


Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}


index_type size()const{return sizeValue;}

index_type dimension(index_type num){return dim[num];}
index_type dimension(index_type num)const{return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){cudaFree(containerMemory);}
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>& inContainer){
rankValue=inContainer.rankValue;
sizeValue=inContainer.sizeValue;
dim[0]=dim0=inContainer.dim0;
dim[1]=dim1=inContainer.dim1;
dim[2]=dim2=inContainer.dim2;
dim[3]=dim3=inContainer.dim3;
dim[4]=dim4=inContainer.dim4;
dim[5]=dim5=inContainer.dim5;
dim[6]=dim6=inContainer.dim6;
dim[7]=dim7=inContainer.dim7;
containerMemory=inContainer.containerMemory;
count_=inContainer.count_;
count_=count_+1;
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}
containerMemory=InContainer.ptr_on_device();

}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutRight,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}


template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>{
index_type dim0=0;
index_type dim1=0;
index_type dim2=0;
index_type dim3=0;
index_type dim4=0;
index_type dim5=0;
index_type dim6=0;
index_type dim7=0;
index_type dim[8]={0};
bool intrepidManaged=false;
Scalar* containerMemory;
unsigned int count_=1;
index_type rankValue=0;
index_type sizeValue=0;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>& InContainer);
FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();


Scalar& operator() (const index_type i0);

Scalar& operator() (const index_type i0, const index_type i1);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 );

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6);

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7);

Scalar& operator() (const index_type i0)const;

Scalar& operator() (const index_type i0, const index_type i1)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 , const index_type i4 )const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6)const;

Scalar& operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4, const index_type i5,
                          const index_type i6, const index_type i7)const;
index_type rank(){return rankValue;}

index_type size(){return sizeValue;}

index_type dimension(index_type num){return dim[num];}

index_type dimension(index_type num)const{return dim[num];}

index_type dimension_0(){return dim0;}
index_type dimension_1(){return dim1;}
index_type dimension_2(){return dim2;}
index_type dimension_3(){return dim3;}
index_type dimension_4(){return dim4;}
index_type dimension_5(){return dim5;}
index_type dimension_6(){return dim6;}
index_type dimension_7(){return dim7;}

index_type dimension_0()const{return dim0;}
index_type dimension_1()const{return dim1;}
index_type dimension_2()const{return dim2;}
index_type dimension_3()const{return dim3;}
index_type dimension_4()const{return dim4;}
index_type dimension_5()const{return dim5;}
index_type dimension_6()const{return dim6;}
index_type dimension_7()const{return dim7;}

};

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intrepidManaged){cudaFree(containerMemory);}
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));

}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}


template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
intrepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
cudaMalloc(&containerMemory,sizeValue*sizeof(Scalar));
}



template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::Rank;
switch(rankValue){
case 1:
sizeValue=dim0;
break;

case 2:
sizeValue=dim0*dim1;
break;

case 3:
sizeValue=dim0*dim1*dim2;
break;

case 4:
sizeValue=dim0*dim1*dim2*dim3;
break;

case 5:
sizeValue=dim0*dim1*dim2*dim3*dim4;
break;

case 6:
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
break;

case 7: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
break;

case 8: 
sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
break;

}

containerMemory=InContainer.ptr_on_device();
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Cuda>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}


#endif
}
#endif
template <class AT>
class randomFillKernel {

 public:
 // typedef typename AT::execution_space execution_space;
 // typedef typename AT::execution_space ExecutionSpace;

//typedef typename Kokkos::OpenMP ExecutionSpace;
//typedef typename Kokkos::OpenMP execution_space; 
 typedef double value_type;
  AT A_, B_, C_;

  randomFillKernel ( AT& A,
               AT& B,
               AT& C) : 
             A_ (A) 
            ,B_ (B)
            ,C_ (C){}

 // typedef Kokkos::View <double * >:: size_type size_type;
  KOKKOS_INLINE_FUNCTION void
  operator () (const index_type i) const
  { // max -plus semiring equivalent of "plus"
    const index_type dim1=A_.dimension(1);
    const index_type dim2=A_.dimension(2);
//  std::default_random_engine generator;
//    std::uniform_real_distribution<double> distribution(0,100);
    
   for (index_type j=0; j< dim1; j++){
    for (index_type k=0; k< dim2; k++){
     A_(i,j,k)=1;
     B_(i,j,k)=2;
     C_(i,j,k)=3;
    }
   }
  }
};

template <class AT>
class mulAddKernel {

 public:
 // typedef typename AT::execution_space execution_space;
 // typedef typename AT::execution_space ExecutionSpace;

// typedef typename Kokkos::OpenMP ExecutionSpace;
//typedef typename Kokkos::OpenMP execution_space; 
 typedef double value_type;
  AT A_, B_, C_;

  mulAddKernel ( AT& A,
              const AT& B,
              const AT& C) :
             A_ (A)
            ,B_ (B)
            ,C_ (C){}

 // typedef Kokkos::View <double * >:: size_type size_type;
  KOKKOS_INLINE_FUNCTION void
  operator () (const index_type i)const
  { // max -plus semiring equivalent of "plus"
   const index_type dim1=A_.dimension(1);
   const index_type dim2=A_.dimension(2);
   for (index_type j=0; j< dim1; j++){
    for (index_type k=0; k< dim2; k++){
     C_(i,j,k)+=B_(i,j,k)*A_(i,j,k);
    }
   }
  }
};
template <class AT>
class TempKernel {

 public:
 // typedef typename AT::execution_space execution_space;
 // typedef typename AT::execution_space ExecutionSpace;

// typedef typename Kokkos::OpenMP ExecutionSpace;
//typedef typename Kokkos::OpenMP execution_space; 
 typedef double value_type;
  AT A_, B_, C_;

  TempKernel ( AT& A,
               AT& B,
               AT& C) :
             A_ (A)
            ,B_ (B)
            ,C_ (C){}

 // typedef Kokkos::View <double * >:: size_type size_type;
  KOKKOS_INLINE_FUNCTION void
  operator () (const index_type i)const
  { // max -plus semiring equivalent of "plus"
   const index_type dim1=A_.dimension(1);
   const index_type dim2=A_.dimension(2);
  //  for (index_type j=0; j< dim1; j++){
  //  for (index_type k=0; k< dim2; k++){
     C_(i,9,9)+=B_(i,9,9)*A_(i,9,9);
  //  }
  // }
  }
};
/*
template <class AT>
__global__ void testcudakernel(AT& A,AT& B, AT& C){

   const index_type dim1=A_.dimension(1);
   const index_type dim2=A_.dimension(2);
    for (index_type j=0; j< dim1; j++){
        for (index_type k=0; k< dim2; k++){
           C_(5,9,9)+=B_(5,9,9)*A_(5,9,9);
              }
           }
  

}
*/

int main(){
index_type problemSize=100;
Kokkos::initialize();
#ifdef KOKKOS_HAVE_SERIAL
{//Test 1 Begin Scope - Kokkos Views Serial vs FieldContainer Kokkos cstd loop

std::cout <<"Test 1 Begin Scope - Kokkos Views Serial vs FieldContainer Kokkos cstd loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Serial>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Serial>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Serial>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


srand(time(NULL));
Kokkos::Impl::Timer RandomViewFillTimer;
for(index_type i=0;i<kokkosviewRightSerialA.dimension_0();i++){
   for(index_type j=0;j<kokkosviewRightSerialA.dimension_1();j++){
      for(index_type k=0;k<kokkosviewRightSerialA.dimension_2();k++){
     kokkosviewRightSerialA(i,j,k)=rand();
     kokkosviewRightSerialB(i,j,k)=rand();
     kokkosviewRightSerialC(i,j,k)=rand();
        }
    }
}
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;
for(index_type i=0;i<fieldcontainerKokkosRightSerialA.dimension_0();i++){
   for(index_type j=0;j<fieldcontainerKokkosRightSerialA.dimension_1();j++){
      for(index_type k=0;k<fieldcontainerKokkosRightSerialA.dimension_2();k++){
     fieldcontainerKokkosRightSerialA(i,j,k)=rand();
     fieldcontainerKokkosRightSerialB(i,j,k)=rand();
     fieldcontainerKokkosRightSerialC(i,j,k)=rand();
        }
    }
}
std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;
Kokkos::fence();
Kokkos::Impl::Timer KokkosForTimer;
for(index_type i=0;i<kokkosviewRightSerialA.dimension_0();i++){
   for(index_type j=0;j<kokkosviewRightSerialA.dimension_1();j++){
      for(index_type k=0;k<kokkosviewRightSerialA.dimension_2();k++){
kokkosviewRightSerialC(i,j,k)+=kokkosviewRightSerialB(i,j,k)*kokkosviewRightSerialA(i,j,k);
       }
    }
}
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
for(index_type i=0;i<fieldcontainerKokkosRightSerialA.dimension_0();i++){
   for(index_type j=0;j<fieldcontainerKokkosRightSerialA.dimension_1();j++){
      for(index_type k=0;k<fieldcontainerKokkosRightSerialA.dimension_2();k++){
fieldcontainerKokkosRightSerialC(i,j,k)+=fieldcontainerKokkosRightSerialB(i,j,k)*fieldcontainerKokkosRightSerialA(i,j,k);
       }
    }
}
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;

}//Test 1 End Scope
#endif

#ifdef KOKKOS_HAVE_OPENMP
{//Test 2 Begin Scope - Kokkos Views OpenMP vs FieldContainer Kokkos cstd loop

std::cout <<"Test 2 Begin Scope - Kokkos Views OpenMP vs FieldContainer Kokkos cstd loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightOpenMPA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightOpenMPB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightOpenMPC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightOpenMPA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightOpenMPB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightOpenMPC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


srand(time(NULL));
Kokkos::Impl::Timer RandomViewFillTimer;
for(index_type i=0;i<kokkosviewRightOpenMPA.dimension_0();i++){
   for(index_type j=0;j<kokkosviewRightOpenMPA.dimension_1();j++){
      for(index_type k=0;k<kokkosviewRightOpenMPA.dimension_2();k++){
     kokkosviewRightOpenMPA(i,j,k)=rand();
     kokkosviewRightOpenMPB(i,j,k)=rand();
     kokkosviewRightOpenMPC(i,j,k)=rand();
        }
    }
}
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;
for(index_type i=0;i<fieldcontainerKokkosRightOpenMPA.dimension_0();i++){
   for(index_type j=0;j<fieldcontainerKokkosRightOpenMPA.dimension_1();j++){
      for(index_type k=0;k<fieldcontainerKokkosRightOpenMPA.dimension_2();k++){
     fieldcontainerKokkosRightOpenMPA(i,j,k)=rand();
     fieldcontainerKokkosRightOpenMPB(i,j,k)=rand();
     fieldcontainerKokkosRightOpenMPC(i,j,k)=rand();
        }
    }
}
std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;
Kokkos::fence();
Kokkos::Impl::Timer KokkosForTimer;
for(index_type i=0;i<kokkosviewRightOpenMPA.dimension_0();i++){
   for(index_type j=0;j<kokkosviewRightOpenMPA.dimension_1();j++){
      for(index_type k=0;k<kokkosviewRightOpenMPA.dimension_2();k++){
kokkosviewRightOpenMPC(i,j,k)+=kokkosviewRightOpenMPB(i,j,k)*kokkosviewRightOpenMPA(i,j,k);
       }
    }
}
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
for(index_type i=0;i<fieldcontainerKokkosRightOpenMPA.dimension_0();i++){
   for(index_type j=0;j<fieldcontainerKokkosRightOpenMPA.dimension_1();j++){
      for(index_type k=0;k<fieldcontainerKokkosRightOpenMPA.dimension_2();k++){
fieldcontainerKokkosRightOpenMPC(i,j,k)+=fieldcontainerKokkosRightOpenMPB(i,j,k)*fieldcontainerKokkosRightOpenMPA(i,j,k);
       }
    }
}
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;

}//Test 2 End Scope
#endif

#ifdef KOKKOS_HAVE_OPENMP
{//Test 3 Begin Scope - Kokkos Views Serial vs FieldContainer Kokkos - Kokkos loop

std::cout <<"Test 3 Begin Scope - Kokkos Views Serial vs FieldContainer Kokkos - Kokkos loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double>fieldcontainerKokkosRightSerialC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


Kokkos::Impl::Timer RandomViewFillTimer;
Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
Kokkos::fill_random(kokkosviewRightSerialA,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialB,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialC,rand_pool64,100);
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;

Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  randomFillKernel< Intrepid2::FieldContainer_Kokkos<double> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );


std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;
Kokkos::fence();


Kokkos::Impl::Timer KokkosForTimer;
Kokkos::parallel_for (kokkosviewRightSerialA.dimension(0),  mulAddKernel< Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::OpenMP> > (kokkosviewRightSerialA, kokkosviewRightSerialB, kokkosviewRightSerialC) );
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  mulAddKernel< Intrepid2::FieldContainer_Kokkos<double> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;

}//Test 3 End Scope
#endif

#ifdef KOKKOS_HAVE_CUDA
{//Test 4 Begin Scope - Kokkos Views CUDA vs FieldContainer Kokkos - Kokkos Managed - Kokkos loop

std::cout <<"Test 4 Begin Scope - Kokkos Views CUDA vs FieldContainer Kokkos - Kokkos Managed - Kokkos loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomViewFillTimer;
Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
Kokkos::fill_random(kokkosviewRightSerialA,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialB,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialC,rand_pool64,100);
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double,double***,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialA(kokkosviewRightSerialA);
Intrepid2::FieldContainer_Kokkos<double,double***,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialB(kokkosviewRightSerialB);
Intrepid2::FieldContainer_Kokkos<double,double***,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialC(kokkosviewRightSerialC);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


Kokkos::Impl::Timer KokkosForTimer;
Kokkos::parallel_for (kokkosviewRightSerialA.dimension(0),  mulAddKernel< Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> > (kokkosviewRightSerialA, kokkosviewRightSerialB, kokkosviewRightSerialC) );
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  mulAddKernel< Intrepid2::FieldContainer_Kokkos<double,double***,Kokkos::LayoutRight,Kokkos::Cuda> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;


}//Test 4 End Scope
#endif

#ifdef KOKKOS_HAVE_CUDA
{//Test 5 Begin Scope - Kokkos Views Cuda vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop

std::cout <<"Test 5 Begin Scope - Kokkos Views Cuda vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda>fieldcontainerKokkosRightSerialC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


Kokkos::Impl::Timer RandomViewFillTimer;
Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
Kokkos::fill_random(kokkosviewRightSerialA,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialB,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialC,rand_pool64,100);
Kokkos::fence();
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  randomFillKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();

std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;

Kokkos::fence();
Kokkos::Impl::Timer KokkosForTimer;
Kokkos::parallel_for (kokkosviewRightSerialA.dimension(0),  mulAddKernel< Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> > (kokkosviewRightSerialA, kokkosviewRightSerialB, kokkosviewRightSerialC) );
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  mulAddKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;

}//Test 5 End Scope
#endif
#ifdef KOKKOS_HAVE_PTHREAD
{//Test 6 Begin Scope - Kokkos Views Threads vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop
std::cout <<"Test 6 Begin Scope - Kokkos Views Threads vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Threads>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Threads>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Threads>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Threads>fieldcontainerKokkosRightSerialA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Threads>fieldcontainerKokkosRightSerialB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Threads>fieldcontainerKokkosRightSerialC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


Kokkos::Impl::Timer RandomViewFillTimer;
Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
Kokkos::fill_random(kokkosviewRightSerialA,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialB,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialC,rand_pool64,100);
Kokkos::fence();
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  randomFillKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Threads> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();

std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;

Kokkos::fence();
Kokkos::Impl::Timer KokkosForTimer;
Kokkos::parallel_for (kokkosviewRightSerialA.dimension(0),  mulAddKernel< Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Threads> > (kokkosviewRightSerialA, kokkosviewRightSerialB, kokkosviewRightSerialC) );
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  mulAddKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Threads> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;
}

#endif

#ifdef KOKKOS_HAVE_CUDA
{//Test 7 Begin Scope - Kokkos Views Cuda vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop

std::cout <<"Test 7 Begin Scope - Kokkos Views Cuda vs FieldContainer Kokkos - Intrepid Managed - Kokkos loop"<<std::endl;
Kokkos::Impl::Timer InitializeViewsTimer;
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialA("A",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialB("B",problemSize,problemSize,problemSize);
Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>kokkosviewRightSerialC("C",problemSize,problemSize,problemSize);
std::cout <<"Initialize Views Time: "<<InitializeViewsTimer.seconds()<<std::endl;

Kokkos::Impl::Timer InitializeFieldContainersTimer;

Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>::execution_space>fieldcontainerKokkosRightSerialA(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>::execution_space>fieldcontainerKokkosRightSerialB(problemSize,problemSize,problemSize);
Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda>::execution_space>fieldcontainerKokkosRightSerialC(problemSize,problemSize,problemSize);
std::cout <<"Initialize FieldContainers Time: "<<InitializeFieldContainersTimer.seconds()<<std::endl;


Kokkos::Impl::Timer RandomViewFillTimer;
Kokkos::Random_XorShift64_Pool<> rand_pool64(5374857);
Kokkos::fill_random(kokkosviewRightSerialA,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialB,rand_pool64,100);
Kokkos::fill_random(kokkosviewRightSerialC,rand_pool64,100);
Kokkos::fence();
std::cout <<"Random View Fill Time: "<<RandomViewFillTimer.seconds()<<std::endl;

Kokkos::Impl::Timer RandomFieldContainerFillTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  randomFillKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();

std::cout <<"Random FieldContainer Fill Time: "<<RandomFieldContainerFillTimer.seconds()<<std::endl;

Kokkos::fence();
Kokkos::Impl::Timer KokkosForTimer;
Kokkos::parallel_for (kokkosviewRightSerialA.dimension(0),  mulAddKernel< Kokkos::View<double***,Kokkos::LayoutRight,Kokkos::Cuda> > (kokkosviewRightSerialA, kokkosviewRightSerialB, kokkosviewRightSerialC) );
Kokkos::fence();
std::cout <<"Kokkos View Time: "<<KokkosForTimer.seconds()<<std::endl;


Kokkos::fence();
Kokkos::Impl::Timer FieldContainerKForTimer;
Kokkos::parallel_for (fieldcontainerKokkosRightSerialA.dimension(0),  mulAddKernel< Intrepid2::FieldContainer_Kokkos<double,void,Kokkos::LayoutRight,Kokkos::Cuda> > (fieldcontainerKokkosRightSerialA, fieldcontainerKokkosRightSerialB, fieldcontainerKokkosRightSerialC) );
Kokkos::fence();
std::cout <<"Field Container Time: "<<FieldContainerKForTimer.seconds()<<std::endl;

}//Test 5 End Scope
#endif

Kokkos::finalize();
return 0;
}

