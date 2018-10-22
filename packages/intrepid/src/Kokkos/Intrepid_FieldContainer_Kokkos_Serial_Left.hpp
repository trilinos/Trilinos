#ifdef KOKKOS_ENABLE_SERIAL
namespace Intrepid{
template <class Scalar,class ScalarPointer>
class FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>{
size_t dim0=0;
size_t dim1=0;
size_t dim2=0;
size_t dim3=0;
size_t dim4=0;
size_t dim5=0;
size_t dim6=0;
size_t dim7=0;
size_t dim[8]={0};
Scalar* containerMemory;
size_t rankValue=0;
size_t sizeValue=0;
bool intrepidManaged=false;
unsigned int count_=1;
public:
FieldContainer_Kokkos()=delete;
FieldContainer_Kokkos(size_t dim_0);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5,size_t dim_6);
FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5,size_t dim_6,size_t dim_7);
FieldContainer_Kokkos(FieldContainer_Kokkos& inContainer);
FieldContainer_Kokkos(const FieldContainer_Kokkos& inContainer);
~FieldContainer_Kokkos();
FieldContainer_Kokkos(Kokkos::View<ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>& InContainer);

typedef Kokkos::Serial execution_space;

Scalar& operator() (const size_t i0);

Scalar& operator() (const size_t i0, const size_t i1);

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2);

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3 );

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3 , const size_t i4 );

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5);

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5,
                          const size_t i6);

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5,
                          const size_t i6, const size_t i7);

Scalar& operator() (const size_t i0)const;

Scalar& operator() (const size_t i0, const size_t i1)const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2)const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3 )const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3 , const size_t i4 )const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5)const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5,
                          const size_t i6)const;

Scalar& operator() (const size_t i0, const size_t i1, const size_t i2,
                          const size_t i3, const size_t i4, const size_t i5,
                          const size_t i6, const size_t i7)const;
size_t rank(){return rankValue;}

size_t size(){return sizeValue;}

size_t dimension(size_t num){return dim[num];}

size_t dimension_0(){return dim0;}
size_t dimension_1(){return dim1;}
size_t dimension_2(){return dim2;}
size_t dimension_3(){return dim3;}
size_t dimension_4(){return dim4;}
size_t dimension_5(){return dim5;}
size_t dimension_6(){return dim6;}
size_t dimension_7(){return dim7;}

size_t dimension_0()const{return dim0;}
size_t dimension_1()const{return dim1;}
size_t dimension_2()const{return dim2;}
size_t dimension_3()const{return dim3;}
size_t dimension_4()const{return dim4;}
size_t dimension_5()const{return dim5;}
size_t dimension_6()const{return dim6;}
size_t dimension_7()const{return dim7;}

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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0){
count_=1;
dim0=dim[0]=dim_0;
rankValue=1;
intrepidManaged=true;
sizeValue=dim_0;
containerMemory=new Scalar[sizeValue];
}

template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1){
count_=1;
dim0=dim[0]=dim_0;
dim1=dim[1]=dim_1;
rankValue=2;
intrepidManaged=true;
sizeValue=dim_0*dim_1;
containerMemory=new Scalar[sizeValue];
}
template <class Scalar,class ScalarPointer>
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2){
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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3){
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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4){
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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5){
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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5,size_t dim_6){
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
FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::FieldContainer_Kokkos(size_t dim_0,size_t dim_1,size_t dim_2,size_t dim_3,size_t dim_4,size_t dim_5,size_t dim_6,size_t dim_7){
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
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0){
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1){
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2){
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3){
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4){
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5){
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5,const size_t i6){
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5,const size_t i6,const size_t i7){
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0)const{
return containerMemory[i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1)const{
return containerMemory[dim0*i1+i0];
}
template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2)const{
return containerMemory[(dim1*i2+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3)const{
return containerMemory[((dim2*i3+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4)const{
return containerMemory[(((dim3*i4+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5)const{
return containerMemory[((((dim4*i5+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5,const size_t i6)const{
return containerMemory[(((((dim5*i6+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}

template <class Scalar,class ScalarPointer>
inline Scalar& FieldContainer_Kokkos<Scalar,ScalarPointer,Kokkos::LayoutLeft,Kokkos::Serial>::operator() (const size_t i0,const size_t i1,const size_t i2,const size_t i3,const size_t i4,const size_t i5,const size_t i6,const size_t i7)const{
return containerMemory[((((((dim6*i7+i6)*dim5+i5)*dim4+i4)*dim3+i3)*dim2+i2)*dim1+i1)*dim0+i0];
}
}
#endif
