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

#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"
#include <Phalanx_any.hpp>
#include "Phalanx_TypeStrings.hpp"
#include <unordered_map>
#include <map>

#include "Sacado.hpp"
#include "Kokkos_View_Fad.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Kokkos_DynRankView_Fad.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#ifdef PHX_ENABLE_KOKKOS_AMT
#include "Kokkos_TaskScheduler.hpp"
#include <type_traits>
#include <limits>
#endif

namespace phalanx_test {

  template <typename Scalar,typename Layout, typename Device, typename... Props>
  class ComputeRhoFlat {
    Kokkos::View<Scalar**,Layout,Device,Props...> rho_;
    Kokkos::View<Scalar**,Layout,Device,Props...> P_;
    Kokkos::View<Scalar**,Layout,Device,Props...> T_;
    Kokkos::View<Scalar*,Layout,Device,Props...> k_;

  public:
    typedef PHX::Device execution_space;

    ComputeRhoFlat(Kokkos::View<Scalar**,Layout,Device,Props...> &rho,
		   Kokkos::View<Scalar**,Layout,Device,Props...> &P,
		   Kokkos::View<Scalar**,Layout,Device,Props...> &T,
		   Kokkos::View<Scalar*,Layout,Device,Props...>& k)
      : rho_(rho)
      , P_(P)
      , T_(T)
      , k_(k) {}

    KOKKOS_INLINE_FUNCTION
    void operator () (const int i) const
    {
      for (int ip = 0; ip < static_cast<int>(rho_.extent(1)); ++ip) {
	rho_(i,ip) = k_(0) * P_(i,ip) / T_(i,ip);
      }
    }
  };

  template <typename Scalar,typename Layout,typename Device, typename... Props>
  class ComputeRhoHierarchic {
    Kokkos::View<Scalar**,Layout,Device,Props...> rho_;
    Kokkos::View<Scalar**,Layout,Device,Props...> P_;
    Kokkos::View<Scalar**,Layout,Device,Props...> T_;
    Kokkos::View<Scalar*,Layout,Device,Props...> k_;

  public:
    typedef PHX::Device execution_space;

    ComputeRhoHierarchic(Kokkos::View<Scalar**,Layout,Device,Props...> &rho,
			 Kokkos::View<Scalar**,Layout,Device,Props...> &P,
			 Kokkos::View<Scalar**,Layout,Device,Props...> &T,
			 Kokkos::View<Scalar*,Layout,Device,Props...>& k)
      : rho_(rho)
      , P_(P)
      , T_(T)
      , k_(k) {}

    KOKKOS_INLINE_FUNCTION
    void operator () (const typename Kokkos::TeamPolicy<Device>::member_type& thread) const
    {
      const int i = thread.league_rank();
      const int num_qp = rho_.extent(1);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(thread,0,num_qp), [=] (const int& ip) {
	rho_(i,ip) = k_(0) * P_(i,ip) / T_(i,ip);
      });
    }
  };

  TEUCHOS_UNIT_TEST(kokkos, FadViewLayoutComparison)
  {
    const int num_cells = 5000;
    const int num_ip = 8;
    // const int deriv_dim = num_ip * 2 + 1;
    const int deriv_dim = 128;

    using FadType = Sacado::Fad::DFad<double>;
    using DefaultLayout = typename PHX::Device::array_layout;
    // using DevLayout = DefaultLayout; // use preferred layout for device
    // using DevLayout = Kokkos::LayoutLeft;
    // using DevLayout = Kokkos::LayoutRight;
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)

#if defined(KOKKOS_ENABLE_CUDA)
    std::cout << "\n\nKOKKOS_ENABLE_CUDA = true" << std::endl;
    const static int FadStride = 32;
#else
    std::cout << "KOKKOS_ENABLE_CUDA = false" << std::endl;
    const static int FadStride = 1;
#endif

    std::cout << "SACADO_VIEW_CUDA_HIERARCHICAL_DFAD = true" << std::endl;
    using DevLayout = Kokkos::LayoutContiguous<DefaultLayout,FadStride>; // use Sacado continguous (best for cuda)
#else
    std::cout << "SACADO_VIEW_CUDA_HIERARCHICAL_DFAD = false" << std::endl;
    using DevLayout = DefaultLayout;
#endif

    std::cout << "DefaultLayout   = " << PHX::typeAsString<DefaultLayout>() << "\n" << std::endl;
    std::cout << "DevLayout   = " << PHX::typeAsString<DevLayout>() << "\n" << std::endl;

    Kokkos::View<FadType**,DevLayout,PHX::Device> rho;
    Kokkos::View<FadType**,DevLayout,PHX::Device> P;
    Kokkos::View<FadType**,DevLayout,PHX::Device> T;
    Kokkos::View<FadType*,DevLayout,PHX::Device> k;
    rho = Kokkos::View<FadType**,DevLayout,PHX::Device>("rho",num_cells,num_ip,deriv_dim);
    P = Kokkos::View<FadType**,DevLayout,PHX::Device>("P",num_cells,num_ip,deriv_dim);
    T = Kokkos::View<FadType**,DevLayout,PHX::Device>("T",num_cells,num_ip,deriv_dim);
    k = Kokkos::View<FadType*,DevLayout,PHX::Device>("k",1,deriv_dim);

    Kokkos::View<FadType**,DevLayout,PHX::Device>::HostMirror host_P;
    Kokkos::View<FadType**,DevLayout,PHX::Device>::HostMirror host_T;
    Kokkos::View<FadType*,DevLayout,PHX::Device>::HostMirror host_k;
    host_P = Kokkos::View<FadType**,DevLayout,PHX::Device>::HostMirror("host_P",num_cells,num_ip,deriv_dim);
    host_T = Kokkos::View<FadType**,DevLayout,PHX::Device>::HostMirror("host_T",num_cells,num_ip,deriv_dim);
    host_k = Kokkos::View<FadType*,DevLayout,PHX::Device>::HostMirror("host_k",1,deriv_dim);


    for (int i=0; i< num_cells; i++){
       for (int j=0; j< num_ip; j++){
         host_P(i,j)=2.0;
	 host_P(i,j).fastAccessDx(j) = 1.0;
         host_T(i,j)=4.0;
	 host_T(i,j).fastAccessDx(num_ip+j) = 1.0;
      }
    }

    host_k(0) = 2.0;
    host_k(0).fastAccessDx(2*num_ip) = 1.0;  // last deriv component is for k

    Kokkos::deep_copy(P, host_P);
    Kokkos::deep_copy(T, host_T);
    Kokkos::deep_copy(k, host_k);
    PHX::Device::fence();

    const int num_samples = 10;

    Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewTimer("Jacobian Flat");
    PHX::Device::fence();
    for (int i=0; i < num_samples; ++i) {
      Teuchos::TimeMonitor tm(*timer);
      Kokkos::parallel_for(num_cells, ComputeRhoFlat<FadType,DevLayout,PHX::Device>(rho, P, T, k));
      PHX::Device::fence();
    }

    timer = Teuchos::TimeMonitor::getNewTimer("Jacobian Hierarchic (AUTO())");
    PHX::Device::fence();
    for (int i=0; i < num_samples; ++i) {
      Teuchos::TimeMonitor tm(*timer);
      Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(num_cells,Kokkos::AUTO()),
			   ComputeRhoHierarchic<FadType,DevLayout,PHX::Device>(rho, P, T, k));
      PHX::Device::fence();
    }

    timer = Teuchos::TimeMonitor::getNewTimer("Jacobian Hierarchic (team=AUTO(),warp=32)");
    PHX::Device::fence();
    for (int i=0; i < num_samples; ++i) {
      Teuchos::TimeMonitor tm(*timer);
      Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(num_cells,Kokkos::AUTO(),32),
			   ComputeRhoHierarchic<FadType,DevLayout,PHX::Device>(rho, P, T, k));
      PHX::Device::fence();
    }

    // ****************************
    // Now redo using the PHX::View
    // ****************************
    {
      PHX::View<FadType**> phx_rho("phx_rho",num_cells,num_ip,deriv_dim);
      PHX::View<FadType**> phx_P("phx_P",num_cells,num_ip,deriv_dim);
      PHX::View<FadType**> phx_T("phx_T",num_cells,num_ip,deriv_dim);
      PHX::View<FadType*> phx_k("phx_k",1,deriv_dim);

      auto phx_host_P = PHX::View<FadType**>::HostMirror("phx_host_P",num_cells,num_ip,deriv_dim);
      auto phx_host_T = PHX::View<FadType**>::HostMirror("phx_host_T",num_cells,num_ip,deriv_dim);
      auto phx_host_k = PHX::View<FadType*>::HostMirror("phx_host_k",1,deriv_dim);

      for (int i=0; i< num_cells; i++){
	for (int j=0; j< num_ip; j++){
	  phx_host_P(i,j).val()=2.0;
	  phx_host_P(i,j).fastAccessDx(j) = 1.0;
	  phx_host_T(i,j).val()=4.0;
	  phx_host_T(i,j).fastAccessDx(num_ip+j) = 1.0;
	}
      }

      phx_host_k(0).val() = 2.0;
      phx_host_k(0).fastAccessDx(2*num_ip) = 1.0;  // last deriv component is for k

      Kokkos::deep_copy(phx_P, phx_host_P);
      Kokkos::deep_copy(phx_T, phx_host_T);
      Kokkos::deep_copy(phx_k, phx_host_k);
      PHX::Device::fence();
      timer = Teuchos::TimeMonitor::getNewTimer("Jacobian Hierarchic <PHX::View> (team=AUTO(),warp=32)");
      PHX::Device::fence();
      for (int i=0; i < num_samples; ++i) {
	Teuchos::TimeMonitor tm(*timer);
	Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(num_cells,Kokkos::AUTO(),32),
			     ComputeRhoHierarchic<FadType,typename PHX::DevLayout<FadType>::type,PHX::Device>(phx_rho, phx_P, phx_T, phx_k));
	PHX::Device::fence();
      }
    }


    std::cout << std::endl;
    Teuchos::TimeMonitor::summarize();
    std::cout << std::endl;
  }

}
