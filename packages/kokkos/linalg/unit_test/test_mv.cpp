#include <cstdio>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <limits>
#include <limits.h>
#include <cmath>

#include <KokkosArray_Host.hpp>
#include <KokkosArray_Cuda.hpp>
#include "KokkosArray_MultiVector.hpp"
#ifndef DEVICE
#define DEVICE 1
#endif
#if DEVICE==1
typedef KokkosArray::Host device_type;
#define KokkosArrayHost(a) a
#define KokkosArrayCUDA(a)
#else
typedef KokkosArray::Cuda device_type;
#define KokkosArrayHost(a)
#define KokkosArrayCUDA(a) a
#endif

#define EPSILON 1e-4


template<typename SCALAR, typename ORDINAL>
int test_mv_dot(ORDINAL size, int numVecs)
{
  typedef typename KokkosArray::MultiVectorDynamic<SCALAR,device_type>::type mv_type;
  typedef typename mv_type::HostMirror h_mv_type;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
  typedef typename vector_type::HostMirror h_vector_type;
  typedef typename mv_type::size_type            size_type;

  mv_type x("X",size,numVecs);
  mv_type y("Y",size,numVecs);
  mv_type r("R",size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = KokkosArray::create_mirror_view(x);
  h_mv_type h_y = KokkosArray::create_mirror_view(y);
  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
  h_vector_type h_a = KokkosArray::create_mirror_view(a);
  h2_vector_type h_b("h2",numVecs);

  srand(17231);
  for(ORDINAL k=0;k<numVecs;k++){
    h_b(k) = 0;
	h_a(k) = 0;
	for(ORDINAL i=0; i<size;i++) {

	  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
	  h_y(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
	  h_b(k)+= h_y(i,k)*h_x(i,k);
	}
  }

  KokkosArray::deep_copy(x,h_x);
  KokkosArray::deep_copy(y,h_y);
  KokkosArray::deep_copy(a,h_a);
  KokkosArray::MV_Dot(a,x,y);
  device_type::fence();

  KokkosArray::deep_copy(h_a,a);
  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k)-h_b(k))/(h_b(k)+EPSILON));
	if(fabs((h_a(k)-h_b(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
  }
  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_add_r_x_y(ORDINAL size, int numVecs){
	  typedef typename KokkosArray::MultiVectorDynamic<ORDINAL,device_type>::type mv_type;
	  typedef typename mv_type::HostMirror h_mv_type;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
	  typedef typename vector_type::HostMirror h_vector_type;
	  typedef typename mv_type::size_type            size_type;

	  mv_type x("X",size,numVecs);
	  mv_type y("Y",size,numVecs);
	  mv_type r("R",size,numVecs);
	  vector_type a("A",numVecs);
	  h_mv_type h_x = KokkosArray::create_mirror_view(x);
	  h_mv_type h_y = KokkosArray::create_mirror_view(y);
	  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
	  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
	  h_vector_type h_a = KokkosArray::create_mirror_view(a);
	  h_vector_type h_b("h_b",numVecs);

	  srand(17231);
	  for(ORDINAL k=0;k<numVecs;k++){
	    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
		for(ORDINAL i=0; i<size;i++) {
		  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_y(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_rh(i,k) = h_y(i,k) + h_x(i,k);
		}
	  }

	  KokkosArray::deep_copy(x,h_x);
	  KokkosArray::deep_copy(y,h_y);
	  KokkosArray::deep_copy(a,h_a);
	  KokkosArray::MV_Add(r,x,y);
	  device_type::fence();

	  KokkosArray::deep_copy(h_rd,r);
	  for(int k=0;k<numVecs;k++){
		h_a(k) = 0;
		h_b(k) = 0;
	    for(int i=0; i<size;i++) {
		  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
		  h_b(k)+= h_rh(i,k)*h_rh(i,k);
		}
	  }

	  double errorsum=0;
	  int errors=0;
	  for(int k=0;k<numVecs;k++)
	  {
	    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
		if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
	  }

	  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_add_r_x_by(ORDINAL size, int numVecs){
	  typedef typename KokkosArray::MultiVectorDynamic<ORDINAL,device_type>::type mv_type;
	  typedef typename mv_type::HostMirror h_mv_type;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
	  typedef typename vector_type::HostMirror h_vector_type;
	  typedef typename mv_type::size_type            size_type;

	  mv_type x("X",size,numVecs);
	  mv_type y("Y",size,numVecs);
	  mv_type r("R",size,numVecs);
	  vector_type a("A",numVecs);
	  h_mv_type h_x = KokkosArray::create_mirror_view(x);
	  h_mv_type h_y = KokkosArray::create_mirror_view(y);
	  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
	  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
	  h_vector_type h_a = KokkosArray::create_mirror_view(a);
	  h_vector_type h_b("h_b",numVecs);

	  srand(17231);
	  for(ORDINAL k=0;k<numVecs;k++){
	    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
		for(ORDINAL i=0; i<size;i++) {
		  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_y(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_rh(i,k) = h_a(k)*h_y(i,k) + h_x(i,k);
		}
	  }

	  KokkosArray::deep_copy(x,h_x);
	  KokkosArray::deep_copy(y,h_y);
	  KokkosArray::deep_copy(a,h_a);
	  KokkosArray::MV_Add(r,x,a,y);
	  device_type::fence();

	  KokkosArray::deep_copy(h_rd,r);
	  for(int k=0;k<numVecs;k++){
		h_a(k) = 0;
		h_b(k) = 0;
	    for(int i=0; i<size;i++) {
		  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
		  h_b(k)+= h_rh(i,k)*h_rh(i,k);
		}
	  }

	  double errorsum=0;
	  int errors=0;
	  for(int k=0;k<numVecs;k++)
	  {
	    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
		if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
	  }

	  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_add_r_ax_by(ORDINAL size, int numVecs){
	  typedef typename KokkosArray::MultiVectorDynamic<ORDINAL,device_type>::type mv_type;
	  typedef typename mv_type::HostMirror h_mv_type;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
	  typedef typename vector_type::HostMirror h_vector_type;
	  typedef typename mv_type::size_type            size_type;

	  mv_type x("X",size,numVecs);
	  mv_type y("Y",size,numVecs);
	  mv_type r("R",size,numVecs);
	  vector_type a("A",numVecs);
	  h_mv_type h_x = KokkosArray::create_mirror_view(x);
	  h_mv_type h_y = KokkosArray::create_mirror_view(y);
	  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
	  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
	  h_vector_type h_a = KokkosArray::create_mirror_view(a);
	  h_vector_type h_b("h_b",numVecs);

	  srand(17231);
	  for(ORDINAL k=0;k<numVecs;k++){
	    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
		for(ORDINAL i=0; i<size;i++) {
		  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_y(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_rh(i,k) = h_a(k)*h_y(i,k) + h_a(k)*h_x(i,k);
		}
	  }

	  KokkosArray::deep_copy(x,h_x);
	  KokkosArray::deep_copy(y,h_y);
	  KokkosArray::deep_copy(a,h_a);
	  KokkosArray::MV_Add(r,a,x,a,y);
	  device_type::fence();

	  KokkosArray::deep_copy(h_rd,r);
	  for(int k=0;k<numVecs;k++){
		h_a(k) = 0;
		h_b(k) = 0;
	    for(int i=0; i<size;i++) {
		  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
		  h_b(k)+= h_rh(i,k)*h_rh(i,k);
		}
	  }

	  double errorsum=0;
	  int errors=0;
	  for(int k=0;k<numVecs;k++)
	  {
	    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
		if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
	  }

	  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_add_r_ax_by(ORDINAL size, int numVecs,int A,int B){
	  typedef typename KokkosArray::MultiVectorDynamic<ORDINAL,device_type>::type mv_type;
	  typedef typename mv_type::HostMirror h_mv_type;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
	  typedef KokkosArray::View<ORDINAL* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
	  typedef typename vector_type::HostMirror h_vector_type;
	  typedef typename mv_type::size_type            size_type;

	  mv_type x("X",size,numVecs);
	  mv_type y("Y",size,numVecs);
	  mv_type r("R",size,numVecs);
	  vector_type a("A",numVecs);
	  h_mv_type h_x = KokkosArray::create_mirror_view(x);
	  h_mv_type h_y = KokkosArray::create_mirror_view(y);
	  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
	  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
	  h_vector_type h_a = KokkosArray::create_mirror_view(a);
	  h_vector_type h_b("h_b",numVecs);

	  srand(17231);
	  for(ORDINAL k=0;k<numVecs;k++){
	    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
		for(ORDINAL i=0; i<size;i++) {
		  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  h_y(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
		  SCALAR aa = A==2?h_a(k):A;
		  SCALAR bb = B==2?h_a(k):B;
		  h_rh(i,k) = bb*h_y(i,k) + aa*h_x(i,k);
		}
	  }

	  KokkosArray::deep_copy(x,h_x);
	  KokkosArray::deep_copy(y,h_y);
	  KokkosArray::deep_copy(a,h_a);
	  KokkosArray::MV_Add(r,a,x,a,y,A,B);
	  device_type::fence();

	  KokkosArray::deep_copy(h_rd,r);
	  for(int k=0;k<numVecs;k++){
		h_a(k) = 0;
		h_b(k) = 0;
	    for(int i=0; i<size;i++) {
		  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
		  h_b(k)+= h_rh(i,k)*h_rh(i,k);
		}
	  }

	  double errorsum=0;
	  int errors=0;
	  for(int k=0;k<numVecs;k++)
	  {
	    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
		if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
	  }

	  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_add(ORDINAL size, int numVecs)
{
  int errors = 0;
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs);
  errors += test_mv_add_r_x_by<SCALAR,ORDINAL>(size,numVecs);
  errors += test_mv_add_r_x_y<SCALAR,ORDINAL>(size,numVecs);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,-1,-1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,-1,1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,-1,2);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,1,-1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,1,1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,1,2);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,2,-1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,2,1);
  errors += test_mv_add_r_ax_by<SCALAR,ORDINAL>(size,numVecs,2,2);
  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_mulscalar_self(ORDINAL size, int numVecs)
{
  typedef typename KokkosArray::MultiVectorDynamic<SCALAR,device_type>::type mv_type;
  typedef typename  mv_type::HostMirror h_mv_type;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
  typedef typename vector_type::HostMirror h_vector_type;
  typedef typename mv_type::size_type            size_type;

  mv_type x("X",(int)size,numVecs);
  mv_type r("R",(int)size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = KokkosArray::create_mirror_view(x);
  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
  h_vector_type h_a = KokkosArray::create_mirror_view(a);
  h_vector_type h_b("h_b",numVecs);

  srand(17231);
  for(int k=0;k<numVecs;k++){
    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
	for(int i=0; i<size;i++) {
	  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
	  h_rh(i,k) = h_a(k)*h_x(i,k);
	}
  }

  KokkosArray::deep_copy(x,h_x);
  KokkosArray::deep_copy(a,h_a);
  KokkosArray::MV_MulScalar(x,a,x);
  device_type::fence();

  KokkosArray::deep_copy(h_rd,x);
  for(int k=0;k<numVecs;k++){
	h_a(k) = 0;
	h_b(k) = 0;
    for(int i=0; i<size;i++) {
	  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
	  h_b(k)+= h_rh(i,k)*h_rh(i,k);
	}
  }

  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
	if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
  }

  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv_mulscalar_diff(ORDINAL size, int numVecs)
{
  typedef typename KokkosArray::MultiVectorDynamic<SCALAR,device_type>::type mv_type;
  typedef typename  mv_type::HostMirror h_mv_type;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,device_type >  vector_type ;
  typedef KokkosArray::View<SCALAR* ,KokkosArray::LayoutLeft,KokkosArray::Host >  h2_vector_type ;
  typedef typename vector_type::HostMirror h_vector_type;
  typedef typename mv_type::size_type            size_type;

  mv_type x("X",(int)size,numVecs);
  mv_type r("R",(int)size,numVecs);
  vector_type a("A",numVecs);
  h_mv_type h_x = KokkosArray::create_mirror_view(x);
  h_mv_type h_rh = KokkosArray::create_mirror_view(r);
  h_mv_type h_rd = KokkosArray::create_mirror_view(r);
  h_vector_type h_a = KokkosArray::create_mirror_view(a);
  h_vector_type h_b("h_b",numVecs);

  srand(17231);
  for(int k=0;k<numVecs;k++){
    h_a(k) = (SCALAR) (1.0*(rand()%40)-20.);
	for(int i=0; i<size;i++) {
	  h_x(i,k) = (SCALAR) (1.0*(rand()%40)-20.);
	  h_rh(i,k) = h_a(k)*h_x(i,k);
	}
  }

  KokkosArray::deep_copy(x,h_x);
  KokkosArray::deep_copy(a,h_a);
  KokkosArray::MV_MulScalar(r,a,x);
  device_type::fence();

  KokkosArray::deep_copy(h_rd,r);
  for(int k=0;k<numVecs;k++){
	h_a(k) = 0;
	h_b(k) = 0;
    for(int i=0; i<size;i++) {
	  h_a(k)+= (h_rh(i,k)-h_rd(i,k))*(h_rh(i,k)-h_rd(i,k));
	  h_b(k)+= h_rh(i,k)*h_rh(i,k);
	}
  }

  double errorsum=0;
  int errors=0;
  for(int k=0;k<numVecs;k++)
  {
    errorsum+=fabs((h_a(k))/(h_b(k)+EPSILON));
	if(fabs((h_a(k))/(h_b(k)+EPSILON))>EPSILON) errors++;
  }

  return errors;
}
template<typename SCALAR, typename ORDINAL>
int test_mv_mulscalar(ORDINAL size, int numVecs)
{
  int errors=0;
  errors += test_mv_mulscalar_self<SCALAR,ORDINAL>(size,numVecs);
  errors += test_mv_mulscalar_diff<SCALAR,ORDINAL>(size,numVecs);
  return errors;
}

template<typename SCALAR, typename ORDINAL>
int test_mv(ORDINAL size, int numVecs, int test){
  if(test==1) return test_mv_dot<SCALAR,ORDINAL>(size,numVecs);
  if(test==2) return test_mv_add<SCALAR,ORDINAL>(size,numVecs);
  if(test==3) return test_mv_mulscalar<SCALAR,ORDINAL>(size,numVecs);
  return 0;
}
template<typename SCALAR, typename ORDINAL>
int test_mv_numVecs(ORDINAL size, int numVecs, int test){
  if(numVecs==-1) {
	int errors = 0;
	for(int nV = 1; nV<20; nV++) {
	  errors+=test_mv<SCALAR,ORDINAL>(size,nV,test);
    }
	return errors;
  }
  return test_mv<SCALAR,ORDINAL>(size,numVecs,test);
}


template<typename ORDINAL>
int test_mv_type(ORDINAL size, int numVecs, int type, int test) {
  if(type==1) return test_mv_numVecs<int,ORDINAL>(size,numVecs,test);
  if(type==2) return test_mv_numVecs<long long int,ORDINAL>(size,numVecs,test);
  if(type==3) return test_mv_numVecs<float,ORDINAL>(size,numVecs,test);
  if(type==4) return test_mv_numVecs<double,ORDINAL>(size,numVecs,test);
  return 0;
}

int test_mv_test(long long int size, int numVecs, int type, int test) {
  int total_errors = 0;
  int mintype = type;
  int maxtype = type;
  if(type<1) {mintype = 1; maxtype = 4;}

  int mintest = test;
  int maxtest = test;
  if(test<1) {mintest = 1; maxtest = 3;}

  type = mintype;
  test = mintest;
  while(test<=maxtest) {
    if(size*numVecs>1<<30) {
      int errors = test_mv_type<long long int>(size,numVecs,type,test);
      printf("Test %s %s with SCALAR=%s and ORDINAL=%s; elements: %li numVecs: %i\n",
	          test==1?"Dot      ":test==2?"Add      ":test==3?"MulScalar":"Invalid  ",
	          errors==0?"PASSED":"FAILED",
   	          type==1?"Int      ":type==2?"Long Long":type==3?"Float    ":type==4?"Double   ":"Invalid  ",
	          "long long int", size, numVecs );
    } else {

      int errors = test_mv_type<int>(size,numVecs,type,test);
      printf("Test %s %s with SCALAR=%s and ORDINAL=%s; elements: %li numVecs: %i\n",
		      test==1?"Dot      ":test==2?"Add      ":test==3?"MulScalar":"Invalid  ",
	          errors==0?"PASSED":"FAILED",
   	          type==1?"Int      ":type==2?"Long Long":type==3?"Float    ":type==4?"Double   ":"Invalid  ",
	          "int          ", size, numVecs );
      total_errors+=errors;
    }
    type++;
    if(type==maxtype+1) {type = mintype; test++;}
  }
  return total_errors;
}

int main(int argc, char **argv)
{
 long long int size = 110503; // a prime number
 int numVecs = -1;
 int threads=1;
 int device = 0;
 int numa=1;
 int test=-1;
 int type=-1;

 for(int i=0;i<argc;i++)
 {
  if((strcmp(argv[i],"-d")==0)) {device=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-s")==0)) {size=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"-v")==0)) {numVecs=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--threads")==0)) {threads=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--numa")==0)) {numa=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--test")==0)) {test=atoi(argv[++i]); continue;}
  if((strcmp(argv[i],"--type")==0)) {type=atoi(argv[++i]); continue;}
 }


 KokkosArrayCUDA(
   KokkosArray::Cuda::SelectDevice select_device(device);
   KokkosArray::Cuda::initialize( select_device );
 )

 if(numa>1 || threads>1)
 {
   KokkosArray::Host::initialize( numa , threads );
 }

 int numVecsList[10] = {1, 2, 3, 4, 5, 8, 11, 15, 16, 17};
 int maxNumVecs = numVecs==-1?17:numVecs;
 int numVecIdx = 0;
 if(numVecs == -1) numVecs = numVecsList[numVecIdx++];

 int total_errors = 0;
 while(numVecs<=maxNumVecs) {
   total_errors += test_mv_test(size,numVecs,type,test);
   if(numVecs<maxNumVecs) numVecs = numVecsList[numVecIdx++];
   else numVecs++;
 }

 if(total_errors == 0)
   printf("Kokkos::MultiVector Test: Passed\n");
 else
   printf("Kokkos::MultiVector Test: Failed\n");


 KokkosArrayCUDA(KokkosArray::Host::finalize();)
 device_type::finalize(  );
}
