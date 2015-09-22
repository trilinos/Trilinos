#ifndef PHX_MD_FIELD_TO_KOKKOS_H
#define PHX_MD_FIELD_TO_KOKKOS_H

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_config.hpp"
#include "Kokkos_View.hpp"

//namepsace PHX {

template <class ScalarT>
struct MDFieldKokkos{
  typedef Kokkos::View<ScalarT********,PHX::Device> t_view8;
  typedef Kokkos::View<ScalarT*******,PHX::Device> t_view7;
  typedef Kokkos::View<ScalarT******,PHX::Device> t_view6;
  typedef Kokkos::View<ScalarT*****,PHX::Device> t_view5;
  typedef Kokkos::View<ScalarT****,PHX::Device> t_view4;
  typedef Kokkos::View<ScalarT***,PHX::Device> t_view3;
  typedef Kokkos::View<ScalarT**,PHX::Device> t_view2;
  typedef Kokkos::View<ScalarT*,PHX::Device> t_view1;
 
  t_view8 view8;
  t_view7 view7;
  t_view6 view6;
  t_view5 view5;
  t_view4 view4;
  t_view3 view3;  
  t_view2 view2;
  t_view1 view1;

  int rank;
  int size;
 
  MDFieldKokkos (){ rank=0;} 
  MDFieldKokkos(char* label, int i0, int i1, int i2, int i3, int i4, int i5, int i6, int i7) {
     view8 = t_view8(label, i0, i1, i2, i3, i4, i5, i6, i7);
     rank=8;
     size=i0*i1*i2*i3*i4*i5*i6*i7;
  } 
  MDFieldKokkos(char* label, int i0, int i1, int i2, int i3, int i4, int i5, int i6) {
     view7 = t_view7(label, i0, i1, i2, i3, i4, i5, i6);
     rank=7;
     size=i0*i1*i2*i3*i4*i5*i6;
  } 
  MDFieldKokkos(char* label, int i0, int i1, int i2, int i3, int i4, int i5) {
     view6 = t_view6(label, i0, i1, i2, i3, i4, i5);
     rank=6;
     size=i0*i1*i2*i3*i4*i5;
  } 
  MDFieldKokkos(char* label, int i0, int i1, int i2, int i3, int i4) {
     view5 = t_view5(label, i0, i1, i2, i3, i4);
     rank=5;
     size=i0*i1*i2*i3*i4;
  } 
  MDFieldKokkos(char* label, int i0, int i1, int i2, int i3) {
     view4 = t_view4(label, i0, i1, i2, i3);
     rank=4;
     size=i0*i1*i2*i3;
  } 
  MDFieldKokkos(char* label, int i0, int i1, int i2) {
     view3 = t_view3(label, i0, i1, i2);
     rank=3;
     size=i0*i1*i2;
  } 
  MDFieldKokkos(char* label, int i0, int i1) {
     view2 = t_view2(label, i0, i1);
     rank=2;
     size=i0*i1;
  } 
  MDFieldKokkos(char* label, int i0) {
     view1 = t_view1(label, i0);
     rank=1;
     size=i0;
  }
//  MDFieldKokkos(const int & r, const t_view8 &v8,  const t_view7 &v7,  const t_view6 &v6,  const t_view4 &v5
//                               const t_view4 &v4,  const t_view3 &v3,  const t_view2 &v2,  const t_view1 &v1):
//                               view8(v8), view7(v7), view6(v6), view5(v5), view4(v4), view3(v3), view2(v2), view1(v1), rank(r) {
//  
//  }
  

     //rank=
     //Do First Touch
/*     if(rank==8)
       Kokkos::parallel_for();
     if(rank==7)
       Kokkos::parallel_for();
     if(rank==6)
       Kokkos::parallel_for();
     if(rank==5)
       Kokkos::parallel_for();
     if(rank==4)
       Kokkos::parallel_for();
     if(rank==3)
       Kokkos::parallel_for();
     if(rank==2)
       Kokkos::parallel_for();
     if(rank==1)
       Kokkos::parallel_for();
 */

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0) { return view1( i0);} 

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1) { return view2( i0, i1);}
  
  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2) { return view3( i0, i1, i2);}

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2, const int& i3) { return view4( i0, i1, i2, i3);}     

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2, const int& i3, const int& i4) { return view5( i0, i1, i2, i3, i4);}

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2, const int& i3, const int& i4, const int& i5) { return view6( i0, i1, i2, i3, i4, i5);}

   KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2, const int& i3, const int& i4, const int& i5, const int& i6) { return view7( i0, i1, i2, i3, i4, i5, i6);}

  KOKKOS_FORCEINLINE_FUNCTION
  ScalarT operator() (const int& i0, const int& i1, const int& i2, const int& i3, const int& i4, const int& i5, const int& i6, const int& i7) {return view8( i0, i1, i2, i3, i4, i5, i6, i7);}

   KOKKOS_FORCEINLINE_FUNCTION
  ScalarT& operator[] (int i0) 
  { 
   if (rank==1)
    return view1(i0);
   if (rank==2)
     return view2(1,0);
   // return view2(i0/view2.dimension_1(),i0%view2.dimension_1());
   if (rank==3)
    return view3(1,0,0);
//    return view3((i0/(view3.dimension_1()*view3.dimension_2()))
//                 (i0/(view3.dimension_2())%dimension_1(),
//                 (i0)%dimension_2());
   if (rank==4)
    return view4(1,0,0,0);
   if (rank==5)
    return view5(1,0,0,0,0);
   if (rank==6)
    return view6(1,0,0,0,0,0);
   if (rank==7)
    return view7(1,0,0,0,0,0,0);
   if (rank==8)
    return view8(1,0,0,0,0,0,0,0);
  } 

 int dimension(int i)
 {
  if (rank==1)
   return view1.dimension_0();

  if (rank==2){
    if (i==0)
      return view2.dimension_0();
    if (i==1)
      return view2.dimension_1();
  }

  if (rank==3){
    if (i==0)
      return view3.dimension_0();
    if (i==1)
      return view3.dimension_1();
    if (i==2)
     return view3.dimension_2();
  }

  if (rank==4){
    if (i==0)
      return view4.dimension_0();
    if (i==1)
      return view4.dimension_1();
    if (i==2)
     return view4.dimension_2();
    if (i==3)
     return view4.dimension_3();
  }

  if (rank==5){
    if (i==0)
      return view5.dimension_0();
    if (i==1)
      return view5.dimension_1();
    if (i==2)
     return view5.dimension_2();
    if (i==3)
     return view5.dimension_3();
    if (i==4)
     return view5.dimension_4(); 
  }

 if (rank==6){
    if (i==0)
      return view6.dimension_0();
    if (i==1)
      return view6.dimension_1();
    if (i==2)
     return view6.dimension_2();
    if (i==3)
     return view6.dimension_3();
    if (i==4)
     return view6.dimension_4();
    if (i==5)
     return view6.dimension_5();
  }
 
  if (rank==7){
    if (i==0)
      return view7.dimension_0();
    if (i==1)
      return view7.dimension_1();
    if (i==2)
     return view7.dimension_2();
    if (i==3)
     return view7.dimension_3();
    if (i==4)
     return view7.dimension_4();
    if (i==5)
     return view7.dimension_5();
    if (i==6)
     return view7.dimension_6();
  }

  if (rank==8){
    if (i==0)
      return view8.dimension_0();
    if (i==1)
      return view8.dimension_1();
    if (i==2)
     return view8.dimension_2();
    if (i==3)
     return view8.dimension_3();
    if (i==4)
     return view8.dimension_4();
    if (i==5)
     return view8.dimension_5();
    if (i==6)
     return view8.dimension_6();
    if (i==7)
     return view8.dimension_7();
  }


 
 }

};
//}
#endif
