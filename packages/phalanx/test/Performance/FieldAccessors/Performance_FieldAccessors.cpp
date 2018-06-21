// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Phalanx_DimTag.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Phalanx_KokkosViewFactory.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_MDField.hpp"
#include "Kokkos_View.hpp"
#include "Shards_Array.hpp"

// From test/Utilities directory
#include "Traits.hpp"

/*! \brief Test to check performance of PHX::MDField vs Kokkos::View

  This test is to show that there is no performance penalty between a
  Kokkos::View and a PHX::MDField.  We wrap a Kokkos::View in the
  MDField.  The compiler should optimize away the wrapper code fro
  optimal performance.  The important comparison is the Compiletime
  MDField and the non-static Kokkos::View<double***>.  At optimization
  of -O3 on gcc 4.8.1, we see virtually no difference in runtime.  We
  also compare against shards arrray to check backwards compatibility
  with original phalanx.  We also look at raw arrays (double*) and
  static Kokkos::Views for an idea of optimal performance.  This test
  is really only valid for Kokkos::Serial nodes such that the amount
  of work done is the same for the objects that do not use the
  "parallel_for" functor.

  NOTE: The actual array sizes used for timings are commented out to
  make the problem smaller.  This is so that this unit test will not
  timeout in debug builds where the executable is much slower.  look
  for the declaration of the num_loops, num_cells, and the static
  Kokkos::View to change the problem size.  These declarations are
  grouped together for convenience.

*/

struct Point : public PHX::DimTag {
  Point(){};
  const char * name() const ;
  static const Point& tag();
};

const char * Point::name() const 
{ static const char n[] = "Point" ; return n ; }
const Point & Point::tag() 
{ static const Point myself ; return myself ; }

struct SPoint : public shards::ArrayDimTag {
  SPoint(){};
  const char * name() const ;
  static const SPoint& tag();
};
const char * SPoint::name() const 
{ static const char n[] = "SPoint" ; return n ; }
const SPoint & SPoint::tag() 
{ static const SPoint myself ; return myself ; }

typedef PHX::index_size_type size_type;

template <typename Scalar,typename Device,typename Array>
class ComputeA {
  Array a_;
  Array b_;
  Array c_;
public:
  typedef PHX::Device execution_space;
  
  ComputeA(Array& a,Array& b,Array& c)
  : a_(a), b_(b), c_(c)
  {}
  
  KOKKOS_INLINE_FUNCTION
  void operator () (const size_type c) const
  {
    const size_type num_ip = a_.extent(1);
    const size_type num_dim = a_.extent(2);
    for (size_type i = 0; i < num_ip; ++i) {
      for (size_type d = 0; d < num_dim; ++d) {
	a_(c,i,d) =  b_(c,i,d) * c_(c,i,d) + b_(c,i,d) + b_(c,i,d) / c_(c,i,d);
      }
    }
  }
};

template <typename Scalar,typename Device,typename Array>
class ComputeARuntime {
  Array a_;
  Array b_;
  Array c_;
public:
  typedef PHX::Device execution_space;
  
  ComputeARuntime(Array& a,Array& b,Array& c)
  : a_(a), b_(b), c_(c)
  {}
  
  KOKKOS_INLINE_FUNCTION
  void operator () (const size_type c) const
  {
    const size_type num_ip = a_.dimension(1);
    const size_type num_dim = a_.dimension(2);
    for (size_type i = 0; i < num_ip; ++i) {
      for (size_type d = 0; d < num_dim; ++d) {
	a_(c,i,d) = b_(c,i,d) * c_(c,i,d) + b_(c,i,d) + b_(c,i,d) / c_(c,i,d);
      }
    }
  }
};

TEUCHOS_UNIT_TEST(performance, ArrayAccessor)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
  RCP<Time> phx_ct_time = TimeMonitor::getNewTimer("MDField Compiletime Rank (no parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> phx_ct_time_pf = TimeMonitor::getNewTimer("MDField Compiletime Rank (with parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> phx_rt_time = TimeMonitor::getNewTimer("MDField Runtime Rank (no parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> phx_rt_time_pf = TimeMonitor::getNewTimer("MDField Runtime Rank (with parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> k_time = TimeMonitor::getNewTimer("KokkosView<double***>(no parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> k_time_pf = TimeMonitor::getNewTimer("KokkosView<double***>(with parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> k_time_static = TimeMonitor::getNewTimer("KokkosView<double[][][]>(no parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> k_time_pf_static = TimeMonitor::getNewTimer("KokkosView<double[][][]>(with parallel_for,device="+PHX::typeAsString<PHX::Device>()+")");
  RCP<Time> s_time = TimeMonitor::getNewTimer("Shards Array Compiletime Rank (no parallel_for)");
  RCP<Time> raw_ptr_time = TimeMonitor::getNewTimer("double* Time");
  
  {    
    TimeMonitor tm_total(*total_time);

    std::cout << std::endl << std::endl
	      << "PHX::Device::size_type = " 
	      << PHX::typeAsString<PHX::Device::size_type>() 
	      << std::endl;

    std::cout << "PHX::index_size_type = " 
	      << PHX::typeAsString<PHX::index_size_type>() 
	      << std::endl;
    
    // For performance testing, build in RELEASE mode and use the
    // numbers in the comments below.  We have small numbers by
    // default for unit testing in debug mode where the code is much
    // slower.
    // const size_type num_loops = 30;
    // const size_type num_cells = 200000;
    // typedef Kokkos::View<double[200000][25][4],PHX::Device> kokkos_field_static;
    const size_type num_loops = 1;
    const size_type num_cells = 100;
    typedef Kokkos::View<double[100][25][4],PHX::Device> kokkos_field_static;
    const size_type num_ip = 25;
    const size_type num_dim = 4;
    const size_type size = num_cells * num_ip * num_dim;

    double* raw_ptr_a = new double[size];
    double* raw_ptr_b = new double[size];
    double* raw_ptr_c = new double[size];
 
    RCP<DataLayout> dl = rcp(new MDALayout<Point,Point,Point>(num_cells,
							      num_ip,
							      num_dim));

    // Compiletime PHX:MDField
    typedef MDField<double,Point,Point,Point> phx_ct_field;
    phx_ct_field phx_ct_a("phx_ct_a", dl);
    phx_ct_field phx_ct_b("phx_ct_a", dl);
    phx_ct_field phx_ct_c("phx_ct_a", dl);
    phx_ct_a.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_ct_a.fieldTag()));
    phx_ct_b.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_ct_b.fieldTag()));
    phx_ct_c.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_ct_c.fieldTag()));

    // Runtime PHX::MDField
    typedef MDField<double> phx_rt_field;
    phx_rt_field phx_rt_a("phx_rt_a", dl);
    phx_rt_field phx_rt_b("phx_rt_b", dl);
    phx_rt_field phx_rt_c("phx_rt_c", dl);
    phx_rt_a.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_rt_a.fieldTag()));
    phx_rt_b.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_rt_b.fieldTag()));
    phx_rt_c.setFieldData(PHX::KokkosViewFactory<double,PHX::Device>::buildView(phx_rt_c.fieldTag()));
    
    // Kokkos View
    typedef Kokkos::View<double***,PHX::Device> kokkos_field;
    kokkos_field kv_a("kv_a",num_cells,num_ip,num_dim);
    kokkos_field kv_b("kv_b",num_cells,num_ip,num_dim);
    kokkos_field kv_c("kv_c",num_cells,num_ip,num_dim);

    // Static Kokkos View
    kokkos_field_static kvs_a("kvs_a");
    kokkos_field_static kvs_b("kvs_b");
    kokkos_field_static kvs_c("kvs_c");
    
    // Compiletime Shards Array
    typedef shards::Array<double,shards::NaturalOrder,SPoint,SPoint,SPoint> shards_field;
    Teuchos::ArrayRCP<double> sa_rcp(size);
    Teuchos::ArrayRCP<double> sb_rcp(size);
    Teuchos::ArrayRCP<double> sc_rcp(size);
    shards_field sa(sa_rcp.get(),num_cells,num_ip,num_dim);
    shards_field sb(sb_rcp.get(),num_cells,num_ip,num_dim);
    shards_field sc(sc_rcp.get(),num_cells,num_ip,num_dim);

    // NOTE: This test assumes a specific layout that may not be
    // optimal for the device.

    for (size_type c=0; c < num_cells; ++c) {
      for (size_type i=0; i < num_ip; ++i) {
	for (size_type d=0; d < num_dim; ++d) {
	  phx_ct_a(c,i,d) = 1.0;
	  phx_ct_b(c,i,d) = 2.0;
	  phx_ct_c(c,i,d) = 3.0;
	  phx_rt_a(c,i,d) = 1.0;
	  phx_rt_b(c,i,d) = 2.0;
	  phx_rt_c(c,i,d) = 3.0;
	  kv_a(c,i,d) = 1.0;
	  kv_b(c,i,d) = 2.0;
	  kv_c(c,i,d) = 3.0;
	  kvs_a(c,i,d) = 1.0;
	  kvs_b(c,i,d) = 2.0;
	  kvs_c(c,i,d) = 3.0;
	  sa(c,i,d) = 1.0;
	  sb(c,i,d) = 2.0;
	  sc(c,i,d) = 3.0;
	  raw_ptr_a[c*(num_ip*num_dim) + i*(num_dim) + d] = 1.0;
	  raw_ptr_b[c*(num_ip*num_dim) + i*(num_dim) + d] = 2.0;
	  raw_ptr_c[c*(num_ip*num_dim) + i*(num_dim) + d] = 3.0;
	}
      }
    }

    // ***********************************************
    // Timings start here
    // ***********************************************
    
    cout << "\nMDField Compiletime Rank (no parallel_for)" << endl;
    {
      TimeMonitor tm(*phx_ct_time);
      for (size_type l=0; l < num_loops; ++l)
	for (size_type c=0; c < num_cells; ++c)
	  for (size_type i=0; i < num_ip; ++i)
	    for (size_type d=0; d < num_dim; ++d)
	      phx_ct_a(c,i,d) = phx_ct_b(c,i,d) * phx_ct_c(c,i,d) + phx_ct_b(c,i,d) + phx_ct_b(c,i,d) / phx_ct_c(c,i,d);
    }

    cout << "MDField Compiletime Rank (with parallel_for)" << endl;
    {
      TimeMonitor tm(*phx_ct_time_pf);
      for (size_type l=0; l < num_loops; ++l) {
	Kokkos::parallel_for(num_cells,ComputeA<double,PHX::Device,phx_ct_field> (phx_ct_a,phx_ct_b,phx_ct_c));
        PHX::Device::fence();
      }
    }
        
    cout << "MDField Runtime Rank (no parallel_for)" << endl;
    {
      TimeMonitor tm(*phx_rt_time);
      for (size_type l=0; l < num_loops; ++l)
	for (size_type c=0; c < num_cells; ++c)
	  for (size_type i=0; i < num_ip; ++i)
	    for (size_type d=0; d < num_dim; ++d)
	      phx_rt_a(c,i,d) = phx_rt_b(c,i,d) * phx_rt_c(c,i,d) + phx_rt_b(c,i,d) + phx_rt_b(c,i,d) / phx_rt_c(c,i,d);
    }

    cout << "MDField Runtime Rank (with parallel_for)" << endl;
    {
      TimeMonitor tm(*phx_rt_time_pf);
      for (size_type l=0; l < num_loops; ++l) {
	Kokkos::parallel_for(num_cells,ComputeARuntime<double,PHX::Device,phx_rt_field> (phx_rt_a,phx_rt_b,phx_rt_c));
        PHX::Device::fence();
      }
    }
    
    cout << "Kokkos View (no parallel_for)" << endl;
    {
      TimeMonitor tm(*k_time);
      for (size_type l=0; l < num_loops; ++l)
	for (size_type c=0; c < num_cells; ++c)
	  for (size_type i=0; i < num_ip; ++i)
	    for (size_type d=0; d < num_dim; ++d)
	      kv_a(c,i,d) = kv_b(c,i,d) * kv_c(c,i,d) + kv_b(c,i,d) + kv_b(c,i,d) / kv_c(c,i,d);
    }
    
    cout << "Kokkos View (parallel_for)" << endl;
    {
      TimeMonitor tm(*k_time_pf);
      for (size_type l=0; l < num_loops; ++l)
	Kokkos::parallel_for(num_cells,ComputeA<double,PHX::Device,kokkos_field>(kv_a,kv_b,kv_c));
      PHX::Device::fence();
    }

    cout << "Static Kokkos View (no parallel_for)" << endl;
    {
      TimeMonitor tm(*k_time_static);
      for (size_type l=0; l <num_loops; ++l) {
	for (size_type c=0; c < num_cells; ++c) {
	  for (size_type i=0; i < num_ip; ++i) {
	    for (size_type d=0; d < num_dim; ++d) {
	      kvs_a(c,i,d) = kvs_b(c,i,d) * kvs_c(c,i,d) + kvs_b(c,i,d) + kvs_b(c,i,d) / kvs_c(c,i,d);
	    }
	  }
	}
      }
    }
    
    cout << "Static Kokkos View (parallel_for)" << endl;
    {
      TimeMonitor tm(*k_time_pf_static);
      for (size_type l=0; l < num_loops; ++l)
	Kokkos::parallel_for(num_cells,ComputeA<double,PHX::Device,kokkos_field_static>(kvs_a,kvs_b,kvs_c));
      PHX::Device::fence();
    }
      
    cout << "Shards Array (no parallel_for)" << endl;
    {
      TimeMonitor tm(*s_time);
      for (size_type l=0; l < num_loops; ++l)
	for (size_type c=0; c < num_cells; ++c)
	  for (size_type i=0; i < num_ip; ++i)
	    for (size_type d=0; d < num_dim; ++d)
	      sa(c,i,d) = sb(c,i,d) * sc(c,i,d) + sb(c,i,d) + sb(c,i,d) / sc(c,i,d);
    }

    cout << "double*" << endl;
    {
      TimeMonitor tm(*raw_ptr_time);
      for (size_type l=0; l < num_loops; ++l)
    	for (size_type c=0; c < num_cells; ++c)
    	  for (size_type i=0; i < num_ip; ++i)
    	    for (size_type d=0; d < num_dim; ++d) {
	      size_type index = c*(num_ip*num_dim) + i*(num_dim) + d;
    	      raw_ptr_a[index] = 
		raw_ptr_b[index] * raw_ptr_c[index] + raw_ptr_b[index] + raw_ptr_b[index] / raw_ptr_c[index];
	    }
    }

    delete [] raw_ptr_a;
    delete [] raw_ptr_b;
    delete [] raw_ptr_c;  
  }

  TimeMonitor::summarize();  
}
