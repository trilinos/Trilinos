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

#ifdef KOKKOS_HAVE_PTHREAD
namespace Intrepid2{
template <class Scalar>
class FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>{
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
bool intepidManaged=true;
index_type count_=1;
public:
FieldContainer_Kokkos(){
count_=1;
dim0=dim[0]=1;
rankValue=1;
intepidManaged=true;
sizeValue=1;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}


FieldContainer_Kokkos(index_type dim_0);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6);
FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7);


void resize(index_type dim_0){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}


dim0=dim[0]=dim_0;
dim1=dim[1]=0;
dim2=dim[2]=0;
dim3=dim[3]=0;
dim4=dim[4]=0;
dim5=dim[5]=0;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=1;
sizeValue=dim_0;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
void resize(index_type dim_0,index_type dim_1){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}

dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=0;
dim3=dim[3]=0;
dim4=dim[4]=0;
dim5=dim[5]=0;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=2;
sizeValue=dim_0*dim_1;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
void resize(index_type dim_0,index_type dim_1,index_type dim_2){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}


dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=0;
dim4=dim[4]=0;
dim5=dim[5]=0;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=3;
sizeValue=dim_0*dim_1*dim_2;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
void resize(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}

dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=0;
dim5=dim[5]=0;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=4;
sizeValue=dim_0*dim_1*dim_2*dim_3;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}


void resize(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}

dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=0;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=5;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
void resize(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}

dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=0;
dim7=dim[7]=0;
rankValue=6;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

void resize(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}

dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=0;
rankValue=7;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;

delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
void resize(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
if(!intepidManaged){
std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
}
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
dim7=dim[7]=dim_7;
rankValue=8;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
delete[] containerMemory;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

void resize(const Teuchos::Array<int>& newDimensions);

template<class Vector>
    void dimensions(Vector& dimensions) const{
 dimensions.resize(rankValue);
 for (index_type i=0; i<rankValue;i++)
   dimensions[i]=dim[i];
}


FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();

typedef Kokkos::Threads execution_space;
  typedef Kokkos::RangePolicy<execution_space,Kokkos::Schedule<Kokkos::Static> > range_policy;

template<class ScalarPoindex_typeer>
FieldContainer_Kokkos(Kokkos::View<ScalarPoindex_typeer,Kokkos::LayoutRight,Kokkos::Threads>& InContainer){
dim0=dim[0]=InContainer.dimension(0);
dim1=dim[1]=InContainer.dimension(1);
dim2=dim[2]=InContainer.dimension(2);
dim3=dim[3]=InContainer.dimension(3);
dim4=dim[4]=InContainer.dimension(4);
dim5=dim[5]=InContainer.dimension(5);
dim6=dim[6]=InContainer.dimension(6);
dim7=dim[7]=InContainer.dimension(7);
rankValue=Kokkos::View<ScalarPoindex_typeer,Kokkos::LayoutRight,Kokkos::Threads>::Rank;
intepidManaged=false;
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

Scalar& operator[] (const index_type i0);

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

Scalar& operator[] (const index_type i0)const;

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
index_type rank() const {return rankValue;}
index_type size(){return sizeValue;}
index_type size() const {return sizeValue;}

index_type dimension(index_type num){return dim[num];}
index_type dimension(index_type num)const {return dim[num];}

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

void initialize(Scalar initValue){

  Kokkos::parallel_for(range_policy(0, sizeValue),initFieldContKokkos<Scalar>(initValue,containerMemory));
}

void initialize(){
Scalar initValue=Scalar(0.0);
 Kokkos::parallel_for(range_policy(0, sizeValue),initFieldContKokkos<Scalar>(initValue,containerMemory));
}

};

template <class Scalar>
void FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::resize(const Teuchos::Array<int>& newDimensions){
// First handle the trivial case of zero dimensions
  if(!intepidManaged){
   std::cerr <<"Resizing Unmanaged FieldContainer_Kokkos Potential Memory Issues"<<std::endl;
   }
  if( newDimensions.size() == 0){
   dim0=dim[0]=0;
   dim1=dim[1]=0;
   dim2=dim[2]=0;
   dim3=dim[3]=0;
   dim4=dim[4]=0;
   dim5=dim[5]=0;
   dim6=dim[6]=0;
   dim7=dim[7]=0;
   rankValue=1;
   sizeValue=0;
   delete[] containerMemory;
   containerMemory=new Scalar[sizeValue];
  }
  else {

    // Copy first 5 dimensions for faster access
    unsigned int theRank = newDimensions.size();
    switch (theRank) {
      case 1:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=0;
        dim2=dim[2]=0;
        dim3=dim[3]=0;
        dim4=dim[4]=0;
        dim5=dim[5]=0;
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 1 ;
        sizeValue=dim0;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

      case 2:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=0;
        dim3=dim[3]=0;
        dim4=dim[4]=0;
        dim5=dim[5]=0;
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 2 ;
        sizeValue=dim0*dim1;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

      case 3:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=0;
        dim4=dim[4]=0;
        dim5=dim[5]=0;
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 3 ;
        sizeValue=dim0*dim1*dim2;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

      case 4:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=newDimensions[3];
        dim4=dim[4]=0;
        dim5=dim[5]=0;
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 4 ;
        sizeValue=dim0*dim1*dim2*dim3;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

      case 5:
       dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=newDimensions[3];
        dim4=dim[4]=newDimensions[4];
        dim5=dim[5]=0;
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 5 ;
        sizeValue=dim0*dim1*dim2*dim3*dim4;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

      case 6:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=newDimensions[3];
        dim4=dim[4]=newDimensions[4];
        dim5=dim[5]=newDimensions[5];
        dim6=dim[6]=0;
        dim7=dim[7]=0;
        rankValue= 6 ;
        sizeValue=dim0*dim1*dim2*dim3*dim4*dim5;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

       case 7:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=newDimensions[3];
        dim4=dim[4]=newDimensions[4];
        dim5=dim[5]=newDimensions[5];
        dim6=dim[6]=newDimensions[6];
        dim7=dim[7]=0;
        rankValue= 7 ;
        sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;

        case 8:
        dim0=dim[0]=newDimensions[0];
        dim1=dim[1]=newDimensions[1];
        dim2=dim[2]=newDimensions[2];
        dim3=dim[3]=newDimensions[3];
        dim4=dim[4]=newDimensions[4];
        dim5=dim[5]=newDimensions[5];
        dim6=dim[6]=newDimensions[6];
        dim7=dim[7]=newDimensions[7];
        rankValue= 8 ;
        sizeValue=dim0*dim1*dim2*dim3*dim4*dim5*dim6*dim7;
        delete[] containerMemory;
        containerMemory=new Scalar[sizeValue];
        break;


      default:
       std::cerr <<"FieldContainer_Kokkos can't have more than 8 dimentions"<<std::endl;
    }
  }
  //this->initialize();
}


template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::~FieldContainer_Kokkos(){
count_=count_-1;
if(count_==0 && intepidManaged){delete[] containerMemory;}
}
template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>& inContainer){
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
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(const FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>& inContainer){
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
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
rankValue=3;
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
rankValue=4;
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
rankValue=5;
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}

template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
rankValue=6;
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}
template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
dim2=dim[2]=dim_2;
dim3=dim[3]=dim_3;
dim4=dim[4]=dim_4;
dim5=dim[5]=dim_5;
dim6=dim[6]=dim_6;
rankValue=7;
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}


template <class Scalar>
FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::FieldContainer_Kokkos(index_type dim_0,index_type dim_1,index_type dim_2,index_type dim_3,index_type dim_4,index_type dim_5,index_type dim_6,index_type dim_7){
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
intepidManaged=true;
sizeValue=dim_0*dim_1*dim_2*dim_3*dim_4*dim_5*dim_6*dim_7;
containerMemory=new Scalar[sizeValue];
//this->initialize();
}


template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator[] (const index_type i0){
return containerMemory[i0];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0){
return containerMemory[i0];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1){
return containerMemory[dim1*i0+i1];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2){
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3){
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4){
return containerMemory[(((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5){
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6){
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7){
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator[] (const index_type i0)const{
return containerMemory[i0];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0)const{
return containerMemory[i0];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1)const{
return containerMemory[dim1*i0+i1];
}
template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2)const{
return containerMemory[(dim1*i0+i1)*dim2+i2];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3)const{
return containerMemory[((dim1*i0+i1)*dim2+i2)*dim3+i3];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4)const{
return containerMemory[(((dim0*i1+i1)*dim2+i2)*dim3+i3)*dim4+i4];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5)const{
return containerMemory[((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6)const{
return containerMemory[(((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6];
}

template <class Scalar>
inline Scalar& FieldContainer_Kokkos<Scalar,Kokkos::LayoutRight,Kokkos::Threads>::operator() (const index_type i0,const index_type i1,const index_type i2,const index_type i3,const index_type i4,const index_type i5,const index_type i6,const index_type i7)const{
return containerMemory[((((((dim1*i0+i1)*dim2+i2)*dim3+i3)*dim4+i4)*dim5+i5)*dim6+i6)*dim7+i7];
}

}
#endif
