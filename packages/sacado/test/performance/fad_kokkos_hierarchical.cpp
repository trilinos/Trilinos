//#define SACADO_VIEW_CUDA_HIERARCHICAL 1
//#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED 1
#define SACADO_VIEW_CUDA_HIERARCHICAL_DFAD 1
#define SACADO_KOKKOS_USE_MEMORY_POOL 1
#define SACADO_ALIGN_SFAD 1

#include "Sacado.hpp"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Kokkos_Core.hpp"
#include "Kokkos_MemoryPool.hpp"
#include "impl/Kokkos_Timer.hpp"
#include <cstdio>
#include <algorithm>

// Advection kernel.

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
struct AdvectionKernel;

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
KOKKOS_INLINE_FUNCTION
AdvectionKernel<FluxView, WgbView, SrcView, WbsView, ResidualView>
create_advection_kernel(const FluxView& flux, const WgbView& bg,
                        const SrcView& src, const WbsView& bs,
                        const ResidualView& residual,
                        const typename FluxView::non_const_value_type& coeff);

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)

// Version of advection kernel that supports flat and hierarchical parallelism
// when using the hierarchical_dfad approach.  This requires no changes to the
// kernel beyond supporting a team policy.
template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
struct AdvectionKernel {
  typedef typename FluxView::non_const_value_type scalar_type;
  typedef typename FluxView::execution_space execution_space;
  typedef typename Kokkos::TeamPolicy<execution_space>::member_type team_handle;

  const FluxView flux_m_i;
  const WgbView wgb;
  const SrcView src_m_i;
  const WbsView wbs;
  const ResidualView residual_m_i;
  const scalar_type coeff;
  const size_t ncells;
  const int num_basis, num_points, num_dim;

  // VS isn't used in this kernel
  template <unsigned VS> struct HierarchicalFlatTag {};
  template <unsigned VS> struct HierarchicalTeamTag {};

  KOKKOS_INLINE_FUNCTION
  AdvectionKernel(const FluxView& flux, const WgbView& gb,
                  const SrcView& src, const WbsView& bs,
                  const ResidualView& residual, const scalar_type& c) :
    flux_m_i(flux),
    wgb(gb),
    src_m_i(src),
    wbs(bs),
    residual_m_i(residual),
    coeff(c),
    ncells(flux_m_i.extent(0)),
    num_basis(wgb.extent(1)),
    num_points(wgb.extent(2)),
    num_dim((wgb.extent(3)))
  {
  }

  KOKKOS_INLINE_FUNCTION
  size_t num_cells() const { return ncells; }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell, const int basis) const {
    scalar_type value(0),value2(0);
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<num_dim; ++dim)
        value += flux_m_i(cell,qp,dim)*wgb(cell,basis,qp,dim);
      value2 += src_m_i(cell,qp)*wbs(cell,basis,qp);
    }
    residual_m_i(cell,basis) = coeff*(value+value2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell) const {
    for (int basis=0; basis<num_basis; ++basis) {
      (*this)(cell,basis);
    }
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalFlatTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank()*team.team_size() + team.team_rank();
    if (cell < ncells)
      (*this)(cell);
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalTeamTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank();
    const int team_size = team.team_size();
    for (int basis=team.team_rank(); basis<num_basis; basis+=team_size)
      (*this)(cell, basis);
  }
};

#elif defined(SACADO_VIEW_CUDA_HIERARCHICAL)

// Version of advection kernel that supports flat and hierarchical parallelism
// when using the hierarchical approach.  This requires a separate scalar type
// for temporaries, and partitioning of non-temporary scalars.
template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
struct AdvectionKernel {
  typedef typename FluxView::non_const_value_type scalar_type;
  typedef typename Kokkos::ThreadLocalScalarType<FluxView>::type local_scalar_type;
  typedef typename FluxView::execution_space execution_space;
  typedef typename Kokkos::TeamPolicy<execution_space>::member_type team_handle;
  enum { stride = Kokkos::ViewScalarStride<FluxView>::stride };

  const FluxView flux_m_i;
  const WgbView wgb;
  const SrcView src_m_i;
  const WbsView wbs;
  const ResidualView residual_m_i;
  const scalar_type coeff;
  const size_t ncells;
  const int num_basis, num_points, num_dim;

  // VS isn't used in this kernel
  template <unsigned VS> struct HierarchicalFlatTag {};
  template <unsigned VS> struct HierarchicalTeamTag {};

  KOKKOS_INLINE_FUNCTION
  AdvectionKernel(const FluxView& flux, const WgbView& gb,
                  const SrcView& src, const WbsView& bs,
                  const ResidualView& residual, const scalar_type& c) :
    flux_m_i(flux),
    wgb(gb),
    src_m_i(src),
    wbs(bs),
    residual_m_i(residual),
    coeff(c),
    ncells(flux_m_i.extent(0)),
    num_basis(wgb.extent(1)),
    num_points(wgb.extent(2)),
    num_dim((wgb.extent(3)))
  {
  }

  KOKKOS_INLINE_FUNCTION
  size_t num_cells() const { return ncells; }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell, const int basis) const {
    local_scalar_type value(0),value2(0);
    local_scalar_type c = Sacado::partition_scalar<stride>(coeff);
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<num_dim; ++dim)
        value += flux_m_i(cell,qp,dim)*wgb(cell,basis,qp,dim);
      value2 += src_m_i(cell,qp)*wbs(cell,basis,qp);
    }
    residual_m_i(cell,basis) = c*(value+value2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell) const {
    for (int basis=0; basis<num_basis; ++basis) {
      (*this)(cell,basis);
    }
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalFlatTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank()*team.team_size() + team.team_rank();
    if (cell < ncells)
      (*this)(cell);
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalTeamTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank();
    const int team_size = team.team_size();
    for (int basis=team.team_rank(); basis<num_basis; basis+=team_size)
      (*this)(cell, basis);
  }
};

#else

// Version of advection kernel that supports flat and hierarchical parallelism
// when using the partitioning approach.  This requires additional code to
// create thread-local views.
template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
struct AdvectionKernel {
  typedef typename FluxView::non_const_value_type scalar_type;
  typedef typename FluxView::execution_space execution_space;
  typedef typename Kokkos::TeamPolicy<execution_space>::member_type team_handle;

  const FluxView flux_m_i;
  const WgbView wgb;
  const SrcView src_m_i;
  const WbsView wbs;
  const ResidualView residual_m_i;
  const scalar_type coeff;
  const size_t ncells;
  const int num_basis, num_points, num_dim;

  template <unsigned VS> struct HierarchicalFlatTag {};
  template <unsigned VS> struct HierarchicalTeamTag {};
  template <unsigned VS> struct PartitionedTag {};

  KOKKOS_INLINE_FUNCTION
  AdvectionKernel(const FluxView& flux, const WgbView& gb,
                  const SrcView& src, const WbsView& bs,
                  const ResidualView& residual, const scalar_type& c) :
    flux_m_i(flux),
    wgb(gb),
    src_m_i(src),
    wbs(bs),
    residual_m_i(residual),
    coeff(c),
    ncells(flux_m_i.extent(0)),
    num_basis(wgb.extent(1)),
    num_points(wgb.extent(2)),
    num_dim((wgb.extent(3)))
  {
  }

  KOKKOS_INLINE_FUNCTION
  size_t num_cells() const { return ncells; }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell, const int basis) const {
    scalar_type value(0),value2(0);
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<num_dim; ++dim)
        value += flux_m_i(cell,qp,dim)*wgb(cell,basis,qp,dim);
      value2 += src_m_i(cell,qp)*wbs(cell,basis,qp);
    }
    residual_m_i(cell,basis) = coeff*(value+value2);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t cell) const {
    for (int basis=0; basis<num_basis; ++basis) {
      (*this)(cell,basis);
    }
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const PartitionedTag<VS>, const size_t cell, const int basis) const {
    // VS is the "vector" size == blockDim.x for CUDA.  This is also set in
    // the team execution policy, but we don't seem to have a way to access it.
    // We also don't have a way to access the vector lane index from the team
    // handle.
#ifdef __CUDA_ARCH__
    const unsigned k = threadIdx.x;
#else
    const unsigned k = 0;
#endif

    // Partition each view based on Cuda thread (vector) index
    auto flux_part =  Kokkos::partition<VS>(flux_m_i, k, VS);
    auto wgb_part =   Kokkos::partition<VS>(wgb, k, VS);
    auto src_part =   Kokkos::partition<VS>(src_m_i, k, VS);
    auto wbs_part =   Kokkos::partition<VS>(wbs, k, VS);
    auto resid_part = Kokkos::partition<VS>(residual_m_i, k, VS);
    auto coeff_part = Sacado::partition_scalar<VS>(coeff);

    // Now run the kernel with thread-local view's
    auto kernel_part = create_advection_kernel(flux_part, wgb_part, src_part,
                                               wbs_part, resid_part,
                                               coeff_part);
    kernel_part(cell, basis);
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalFlatTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank()*team.team_size() + team.team_rank();
    if (cell < ncells)
      for (int basis=0; basis<num_basis; ++basis)
        (*this)(PartitionedTag<VS>(), cell, basis);
  }

  template <unsigned VS>
  KOKKOS_INLINE_FUNCTION
  void operator() (const HierarchicalTeamTag<VS>, const team_handle& team) const {
    const size_t cell = team.league_rank();
    const int team_size = team.team_size();
    for (int basis=team.team_rank(); basis<num_basis; basis+=team_size)
      (*this)(PartitionedTag<VS>(), cell, basis);
  }
};

#endif

template<typename FluxView, typename WgbView, typename SrcView,
         typename WbsView, typename ResidualView>
KOKKOS_INLINE_FUNCTION
AdvectionKernel<FluxView, WgbView, SrcView, WbsView, ResidualView>
create_advection_kernel(const FluxView& flux, const WgbView& bg,
                        const SrcView& src, const WbsView& bs,
                        const ResidualView& residual,
                        const typename FluxView::non_const_value_type& coeff)
{
  typedef AdvectionKernel<FluxView, WgbView, SrcView, WbsView, ResidualView> kernel_type;
  return kernel_type(flux,bg,src,bs,residual,coeff);
}

template<typename KernelType>
void run_flat(const KernelType& kernel) {
  typedef typename KernelType::execution_space execution_space;
  Kokkos::RangePolicy<execution_space> policy(0,kernel.num_cells());
  Kokkos::parallel_for(policy, kernel);
}

template<typename KernelType>
void run_hierarchical_flat(const KernelType& kernel) {
  typedef typename KernelType::execution_space execution_space;
#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
#else
  const bool is_cuda = false;
#endif
  const unsigned vector_size = is_cuda ? 32 : 1;
  if (is_cuda) {
    const unsigned team_size = 256 / vector_size;
    typedef typename KernelType::template HierarchicalFlatTag<vector_size> tag_type;
    typedef Kokkos::TeamPolicy<execution_space,tag_type> policy_type;
    const size_t range = (kernel.num_cells()+team_size-1)/team_size;
    policy_type policy(range,team_size,vector_size);
    Kokkos::parallel_for(policy, kernel);
  }
  else {
    run_flat(kernel);
  }
}

template<typename KernelType>
void run_hierarchical_team(const KernelType& kernel) {
  typedef typename KernelType::execution_space execution_space;
#if defined (KOKKOS_ENABLE_CUDA)
  const bool is_cuda = std::is_same<execution_space, Kokkos::Cuda>::value;
#else
  const bool is_cuda = false;
#endif
  const unsigned vector_size = is_cuda ? 32 : 1;
  if (is_cuda) {
    const unsigned team_size = 256 / vector_size;
    typedef typename KernelType::template HierarchicalTeamTag<vector_size> tag_type;
    typedef Kokkos::TeamPolicy<execution_space,tag_type> policy_type;
    policy_type policy(kernel.num_cells(),team_size,vector_size);
    Kokkos::parallel_for(policy, kernel);
  }
  else {
    run_flat(kernel);
  }
}

template <typename T> struct FadTypeName;
template <typename T> struct FadTypeName< Sacado::Fad::DFad<T> > {
  static std::string eval() {
    return std::string("dfad")
#if defined(SACADO_KOKKOS_USE_MEMORY_POOL)
      + std::string(", mempool")
#endif
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD_STRIDED)
      + std::string(", strided")
#endif
      ;
  }
};
template <typename T, int N> struct FadTypeName< Sacado::Fad::SFad<T,N> > {
  static std::string eval() {
#if defined(SACADO_ALIGN_SFAD)
    return "sfad, aligned";
#else
    return "sfad";
#endif
  }
};
template <typename T, int N> struct FadTypeName< Sacado::Fad::SLFad<T,N> > {
  static std::string eval() { return "slfad"; }
};

template<typename ExecSpace, int DIM, int N>
struct DrekarTest {

  int ncells;
  int num_basis;
  int num_points;
  int ndim;

  struct MomFluxTag {};
  struct MomFluxTagConst {};
  struct MomFluxTagConstTeam {};

  typedef Kokkos::View<double****[N+1],ExecSpace> t_4DView;
  typedef Kokkos::View<double***[N+1],ExecSpace> t_3DView;
  typedef Kokkos::View<double**[N+1],ExecSpace> t_2DView;

  typedef Kokkos::View<const double***[N+1],ExecSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_3DView_const;

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  typedef Sacado::Fad::DFad<double> FadType; // Must be DFad in this case
#else
  typedef Sacado::Fad::SFad<double,N> FadType;
#endif
  typedef Kokkos::View<FadType****,ExecSpace> t_4DViewFad;
  typedef Kokkos::View<FadType***,ExecSpace> t_3DViewFad;
  typedef Kokkos::View<FadType**,ExecSpace> t_2DViewFad;

  typedef typename ExecSpace::array_layout DefaultLayout;
#if defined(KOKKOS_ENABLE_CUDA)
  static const int FadStride =
    std::is_same< ExecSpace, Kokkos::Cuda >::value ? 32 : 1;
#if defined(SACADO_ALIGN_SFAD)
  static const int Nalign = ((N+FadStride-1)/FadStride)*FadStride;
  typedef typename FadType::template apply_N<Nalign>::type AlignedFadType;
#else
  typedef FadType AlignedFadType;
#endif
#else
  static const int FadStride = 1;
  typedef FadType AlignedFadType;
#endif
  typedef Kokkos::LayoutContiguous<DefaultLayout,FadStride> ContLayout;
  typedef Kokkos::View<AlignedFadType****,ContLayout,ExecSpace> t_4DViewFadCont;
  typedef Kokkos::View<AlignedFadType***,ContLayout,ExecSpace> t_3DViewFadCont;
  typedef Kokkos::View<AlignedFadType**,ContLayout,ExecSpace> t_2DViewFadCont;

  typedef typename Kokkos::TeamPolicy<ExecSpace>::member_type team_handle;

  typedef Kokkos::View<double****[N+1],Kokkos::LayoutRight,ExecSpace> t_4DView_team;
  typedef Kokkos::View<double***[N+1],Kokkos::LayoutRight,ExecSpace> t_3DView_team;
  typedef Kokkos::View<double**[N+1],Kokkos::LayoutRight,ExecSpace> t_2DView_team;
  typedef Kokkos::View<const double***[N+1],Kokkos::LayoutRight,ExecSpace,Kokkos::MemoryTraits<Kokkos::RandomAccess> > t_3DView_const_team;


  typedef Kokkos::View<double[N+1],typename ExecSpace::scratch_memory_space,Kokkos::MemoryTraits<Kokkos::Unmanaged> > t_shared_scalar;

  t_4DViewFad wgb_fad;
  t_3DViewFad flux_m_i_fad,wbs_fad;
  t_2DViewFad src_m_i_fad,residual_m_i_fad;

  t_4DViewFadCont wgb_fad_cont;
  t_3DViewFadCont flux_m_i_fad_cont,wbs_fad_cont;
  t_2DViewFadCont src_m_i_fad_cont,residual_m_i_fad_cont;

  t_4DView wgb;
  t_3DView flux_m_i,wbs;
  t_3DView_const flux_m_i_const;
  t_2DView src_m_i,residual_m_i,residual_m_i_const;

  t_4DView_team twgb;
  t_3DView_team tflux_m_i,twbs;
  t_3DView_const_team tflux_m_i_const;
  t_2DView_team tsrc_m_i,tresidual_m_i;

  AlignedFadType coeff;

  DrekarTest(int ncells_, int num_basis_, int num_points_):
    ncells(ncells_) ,
    num_basis(num_basis_) ,
    num_points(num_points_) ,
    ndim(DIM),
    coeff(N, 0.0)
  {
    wgb_fad = t_4DViewFad("",ncells,num_basis,num_points,ndim,N+1);
    wbs_fad = t_3DViewFad("",ncells,num_basis,num_points,N+1);
    flux_m_i_fad = t_3DViewFad("",ncells,num_points,ndim,N+1);
    src_m_i_fad = t_2DViewFad("",ncells,num_points,N+1);
    residual_m_i_fad = t_2DViewFad("",ncells,num_basis,N+1);
    init_fad(wgb_fad, wbs_fad, flux_m_i_fad, src_m_i_fad, residual_m_i_fad);

    wgb_fad_cont = t_4DViewFadCont("",ncells,num_basis,num_points,ndim,N+1);
    wbs_fad_cont = t_3DViewFadCont("",ncells,num_basis,num_points,N+1);
    flux_m_i_fad_cont = t_3DViewFadCont("",ncells,num_points,ndim,N+1);
    src_m_i_fad_cont = t_2DViewFadCont("",ncells,num_points,N+1);
    residual_m_i_fad_cont = t_2DViewFadCont("",ncells,num_basis,N+1);
    init_fad(wgb_fad_cont, wbs_fad_cont, flux_m_i_fad_cont, src_m_i_fad_cont, residual_m_i_fad_cont);

    wgb = t_4DView("",ncells,num_basis,num_points,ndim);
    wbs = t_3DView("",ncells,num_basis,num_points);
    flux_m_i_const = flux_m_i = t_3DView("",ncells,num_points,ndim);
    src_m_i = t_2DView("",ncells,num_points);
    residual_m_i = t_2DView("",ncells,num_basis);
    init_array(wgb, wbs, flux_m_i, src_m_i, residual_m_i);

    residual_m_i_const = t_2DView("",ncells,num_basis);
    Kokkos::deep_copy( residual_m_i_const, 0.0 );

    twgb = t_4DView_team("",ncells,num_basis,num_points,ndim);
    twbs = t_3DView_team("",ncells,num_basis,num_points);
    tflux_m_i_const = tflux_m_i = t_3DView_team("",ncells,num_points,ndim);
    tsrc_m_i = t_2DView_team("",ncells,num_points);
    tresidual_m_i = t_2DView_team("",ncells,num_basis);
    init_array(twgb, twbs, tflux_m_i, tsrc_m_i, tresidual_m_i);

    for (int i=0; i<N; ++i)
      coeff.fastAccessDx(i) = generate_fad(1,1,1,1,N,0,0,0,0,i);
    coeff.val() = generate_fad(1,1,1,1,N,0,0,0,0,N);
  }

  typename FadType::value_type
  generate_fad( const size_t n0, const size_t n1,
                const size_t n2, const size_t n3, const int fad_size,
                const size_t i0, const size_t i1,
                const size_t i2, const size_t i3,
                const int i_fad)
  {
    typedef typename FadType::value_type scalar;
    const scalar x0 =    10.0 + scalar(n0) / scalar(i0+1);
    const scalar x1 =   100.0 + scalar(n1) / scalar(i1+1);
    const scalar x2 =  1000.0 + scalar(n2) / scalar(i2+1);
    const scalar x3 = 10000.0 + scalar(n3) / scalar(i3+1);
    const scalar x  = x0 + x1 + x2 + x3;
    if (i_fad == fad_size)
      return x;
    const scalar x_fad = 1.0 + scalar(fad_size) / scalar(i_fad+1);
    return x + x_fad;
  }

  template <typename V1, typename V2, typename V3, typename V4, typename V5>
  void init_fad(const V1& v1, const V2& v2, const V3& v3, const V4& v4, const V5& v5)
  {
    // Kokkos::deep_copy(typename V1::array_type(v1), 1.0);
    // Kokkos::deep_copy(typename V2::array_type(v2), 2.0);
    // Kokkos::deep_copy(typename V3::array_type(v3), 3.0);
    // Kokkos::deep_copy(typename V4::array_type(v4), 4.0);

    auto v1_h = Kokkos::create_mirror_view(v1);
    auto v2_h = Kokkos::create_mirror_view(v2);
    auto v3_h = Kokkos::create_mirror_view(v3);
    auto v4_h = Kokkos::create_mirror_view(v4);
    for (int cell=0; cell<ncells; ++cell) {
      for (int basis=0; basis<num_basis; ++basis) {
        for (int qp=0; qp<num_points; ++qp) {
          for (int dim=0; dim<ndim; ++dim) {
            for (int i=0; i<N; ++i)
              v1_h(cell,basis,qp,dim).fastAccessDx(i) =
                generate_fad(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,i);
            v1_h(cell,basis,qp,dim).val() =
              generate_fad(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,N);
          }
          for (int i=0; i<N; ++i)
            v2_h(cell,basis,qp).fastAccessDx(i) =
              generate_fad(ncells,num_basis,num_points,1,N,cell,basis,qp,0,i);
          v2_h(cell,basis,qp).val() =
              generate_fad(ncells,num_basis,num_points,1,N,cell,basis,qp,0,N);
        }
      }
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<ndim; ++dim) {
          for (int i=0; i<N; ++i)
            v3_h(cell,qp,dim).fastAccessDx(i) =
              generate_fad(ncells,1,num_points,ndim,N,cell,0,qp,dim,i);
          v3_h(cell,qp,dim).val() =
              generate_fad(ncells,1,num_points,ndim,N,cell,0,qp,dim,N);
        }
        for (int i=0; i<N; ++i)
          v4_h(cell,qp).fastAccessDx(i) =
            generate_fad(ncells,1,num_points,1,N,cell,0,qp,0,i);
        v4_h(cell,qp).val() =
            generate_fad(ncells,1,num_points,1,N,cell,0,qp,0,N);
      }
    }

    Kokkos::deep_copy( v1, v1_h );
    Kokkos::deep_copy( v2, v2_h );
    Kokkos::deep_copy( v3, v3_h );
    Kokkos::deep_copy( v4, v4_h );

    Kokkos::deep_copy(typename V5::array_type(v5), 0.0);
  }

  template <typename V1, typename V2, typename V3, typename V4, typename V5>
  void init_array(const V1& v1, const V2& v2, const V3& v3, const V4& v4, const V5& v5)
  {
    // Kokkos::deep_copy(typename V1::array_type(v1), 1.0);
    // Kokkos::deep_copy(typename V2::array_type(v2), 2.0);
    // Kokkos::deep_copy(typename V3::array_type(v3), 3.0);
    // Kokkos::deep_copy(typename V4::array_type(v4), 4.0);

    auto v1_h = Kokkos::create_mirror_view(v1);
    auto v2_h = Kokkos::create_mirror_view(v2);
    auto v3_h = Kokkos::create_mirror_view(v3);
    auto v4_h = Kokkos::create_mirror_view(v4);
    for (int cell=0; cell<ncells; ++cell) {
      for (int basis=0; basis<num_basis; ++basis) {
        for (int qp=0; qp<num_points; ++qp) {
          for (int dim=0; dim<ndim; ++dim) {
            for (int i=0; i<N; ++i)
              v1_h(cell,basis,qp,dim,i) =
                generate_fad(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,i);
            v1_h(cell,basis,qp,dim,N) =
              generate_fad(ncells,num_basis,num_points,ndim,N,cell,basis,qp,dim,N);
          }
          for (int i=0; i<N; ++i)
            v2_h(cell,basis,qp,i) =
              generate_fad(ncells,num_basis,num_points,1,N,cell,basis,qp,0,i);
          v2_h(cell,basis,qp,N) =
            generate_fad(ncells,num_basis,num_points,1,N,cell,basis,qp,0,N);
        }
      }
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<ndim; ++dim) {
          for (int i=0; i<N; ++i)
            v3_h(cell,qp,dim,i) =
              generate_fad(ncells,1,num_points,ndim,N,cell,0,qp,dim,i);
          v3_h(cell,qp,dim,N) =
            generate_fad(ncells,1,num_points,ndim,N,cell,0,qp,dim,N);
        }
        for (int i=0; i<N; ++i)
          v4_h(cell,qp,i) =
            generate_fad(ncells,1,num_points,1,N,cell,0,qp,0,i);
        v4_h(cell,qp,N) =
          generate_fad(ncells,1,num_points,1,N,cell,0,qp,0,N);
      }
    }

    Kokkos::deep_copy( v1, v1_h );
    Kokkos::deep_copy( v2, v2_h );
    Kokkos::deep_copy( v3, v3_h );
    Kokkos::deep_copy( v4, v4_h );

    Kokkos::deep_copy(typename V5::array_type(v5), 0.0);
  }

  template <typename View1, typename View2>
  typename std::enable_if< !Kokkos::is_view_fad<View2>::value, bool>::type
  check(const View1& v_gold, const View2& v, const double tol)
  {
    // Copy to host
    typename View1::HostMirror v_gold_h = Kokkos::create_mirror_view(v_gold);
    typename View2::HostMirror v_h      = Kokkos::create_mirror_view(v);
    Kokkos::deep_copy(v_gold_h, v_gold);
    Kokkos::deep_copy(v_h, v);

    typedef typename View1::value_type value_type;

    const size_t n0 = v_gold_h.extent(0);
    const size_t n1 = v_gold_h.extent(1);
    const size_t n2 = v_gold_h.extent(2);

    bool success = true;
    for ( size_t i0 = 0 ; i0 < n0 ; ++i0 ) {
      for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
        for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
          value_type x_gold = v_gold_h(i0,i1,i2);
          value_type x      = v_h(i0,i1,i2);
          if (std::abs(x_gold-x) > tol*std::abs(x_gold)) {
            std::cout << "Comparison failed!  x_gold("
                      << i0 << "," << i1 << "," << i2 << ") = "
                      << x_gold << " , x = " << x
                      << std::endl;
            success = false;
          }
        }
      }
    }

    return success;
  }

  template <typename View1, typename View2>
  typename std::enable_if< Kokkos::is_view_fad<View2>::value, bool>::type
  check(const View1& v_gold, const View2& v, const double tol)
  {
    // Copy to host
    typename View1::HostMirror v_gold_h = Kokkos::create_mirror_view(v_gold);
    typename View2::HostMirror v_h      = Kokkos::create_mirror_view(v);
    Kokkos::deep_copy(v_gold_h, v_gold);
    Kokkos::deep_copy(v_h, v);

    typedef typename View1::value_type value_type;

    const size_t n0 = v_gold_h.extent(0);
    const size_t n1 = v_gold_h.extent(1);
    const size_t n2 = v_gold_h.extent(2);

    bool success = true;
    for ( size_t i0 = 0 ; i0 < n0 ; ++i0 ) {
      for ( size_t i1 = 0 ; i1 < n1 ; ++i1 ) {
        for ( size_t i2 = 0 ; i2 < n2 ; ++i2 ) {
          value_type x_gold = v_gold_h(i0,i1,i2);
          value_type x      = (i2 == n2-1) ? v_h(i0,i1).val() : v_h(i0,i1).dx(i2);
          if (std::abs(x_gold-x) > tol*std::abs(x_gold)) {
            std::cout << "Comparison failed!  x_gold("
                      << i0 << "," << i1 << "," << i2 << ") = "
                      << x_gold << " , x = " << x
                      << std::endl;
            success = false;
          }
        }
      }
    }

    return success;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const MomFluxTag, const std::size_t &cell) const {
    for (int basis=0; basis<num_basis; ++basis) {
      double value[N+1],value2[N+1];
      for (int k=0; k<N+1; ++k) {
        value[k] = 0.0;
        value2[k] = 0.0;
      }
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<DIM; ++dim) {
          const double flux_val = flux_m_i(cell,qp,dim,N);
          const double wgb_val = wgb(cell,basis,qp,dim,N);
          value[N] += flux_val*wgb_val;
          for(int k = 0; k<N;k++)
            value[k] += flux_val*wgb(cell,basis,qp,dim,k)+flux_m_i(cell,qp,dim,k)*wgb_val;
        }
        const double src_val = src_m_i(cell,qp,N);
        const double wbs_val = wbs(cell,basis,qp,N);
        value2[N] += src_val*wbs_val;
        for(int k = 0; k<N;k++)
          value2[k] += src_val*wbs(cell,basis,qp,k)+src_m_i(cell,qp,k)*wbs_val;
      }
      for(int k = 0; k<N; k++)
        residual_m_i(cell,basis,k) =
          coeff.val()*(value[k]+value2[k]) +
          coeff.fastAccessDx(k)*(value[N]+value2[N]);
      residual_m_i(cell,basis,N)= coeff.val()*(value[N]+value2[N]);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const MomFluxTagConst, const std::size_t &cell) const {
    for (int basis=0; basis<num_basis; ++basis) {
      double value[N+1],value2[N+1];
      for (int k=0; k<N+1; ++k) {
        value[k] = 0.0;
        value2[k] = 0.0;
      }
      for (int qp=0; qp<num_points; ++qp) {
        for (int dim=0; dim<DIM; ++dim) {
          const double flux_val = flux_m_i(cell,qp,dim,N);
          const double wgb_val = wgb(cell,basis,qp,dim,N);
          value[N] += flux_val*wgb_val;
          for(int k = 0; k<N;k++)
            value[k] += flux_val*wgb(cell,basis,qp,dim,k)+flux_m_i_const(cell,qp,dim,k)*wgb_val;
        }
        const double src_val = src_m_i(cell,qp,N);
        const double wbs_val = wbs(cell,basis,qp,N);
        value2[N] += src_val*wbs_val;
        for(int k = 0; k<N;k++)
          value2[k] += src_val*wbs(cell,basis,qp,k)+src_m_i(cell,qp,k)*wbs_val;
      }
      for(int k = 0; k<N; k++)
        residual_m_i_const(cell,basis,k) =
          coeff.val()*(value[k]+value2[k]) +
          coeff.fastAccessDx(k)*(value[N]+value2[N]);
      residual_m_i_const(cell,basis,N)= coeff.val()*(value[N]+value2[N]);
    }
  }


  KOKKOS_INLINE_FUNCTION
  void compute_one(const MomFluxTagConstTeam, const team_handle& team, const int &cell, const int& basis,
                   const t_shared_scalar& value, const t_shared_scalar& value2) const {
    for (int qp=0; qp<num_points; ++qp) {
      for (int dim=0; dim<DIM; ++dim) {
        const double flux_val = tflux_m_i(cell,qp,dim,N);
        const double wgb_val = twgb(cell,basis,qp,dim,N);
        Kokkos::single(Kokkos::PerThread(team), [&] () {
            value[N] += flux_val*wgb_val;
          });
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N), [&] (const int& k) {
            value[k] += flux_val*twgb(cell,basis,qp,dim,k)+tflux_m_i_const(cell,qp,dim,k)*wgb_val;
          });
      }
      const double src_val = tsrc_m_i(cell,qp,N);
      const double wbs_val = twbs(cell,basis,qp,N);
      Kokkos::single(Kokkos::PerThread(team), [&] () {
          value2[N] += src_val*wbs_val;
        });
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N), [&] (const int& k) {
          value2[k] += src_val*twbs(cell,basis,qp,k)+tsrc_m_i(cell,qp,k)*wbs_val;
        });
    }
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N), [&] (const int& k) {
        tresidual_m_i(cell,basis,k) =
          coeff.val()*(value[k]+value2[k]) +
          coeff.fastAccessDx(k)*(value[N]+value2[N]);
      });
    Kokkos::single(Kokkos::PerThread(team), [&] () {
          tresidual_m_i(cell,basis,N) = coeff.val()*(value[N]+value2[N]);
        });
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const MomFluxTagConstTeam, const team_handle& team) const {
    t_shared_scalar value1(team.thread_scratch(0));
    t_shared_scalar value2(team.thread_scratch(0));

    const int cell = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_basis), [&] (const int& basis) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,N+1), [&] (const int& k) {
            value1(k) = 0;
            value2(k) = 0;
          });
        compute_one(MomFluxTagConstTeam(),team,cell,basis,value1,value2);
      });
  }

  void compute(const int ntrial, const bool do_check) {
    auto kernel_flat =
      create_advection_kernel(flux_m_i_fad, wgb_fad, src_m_i_fad, wbs_fad, residual_m_i_fad, coeff);
    Kokkos::Impl::Timer timer;
    for (int i=0; i<ntrial; ++i) {
      run_flat(kernel_flat);
    }
    Kokkos::fence();
    double time_fad = timer.seconds() / ntrial / ncells;

    auto kernel_team =
      create_advection_kernel(flux_m_i_fad_cont, wgb_fad_cont, src_m_i_fad_cont, wbs_fad_cont, residual_m_i_fad_cont, coeff);
    timer.reset();
    for (int i=0; i<ntrial; ++i) {
      run_hierarchical_team(kernel_team);
    }
    Kokkos::fence();
    double time_fad_cont = timer.seconds() / ntrial / ncells;

    timer.reset();
    for (int i=0; i<ntrial; ++i)
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace,MomFluxTag>(0,ncells), *this);
    Kokkos::fence();
    double time = timer.seconds() / ntrial / ncells;

    timer.reset();
    for (int i=0; i<ntrial; ++i)
      Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace,MomFluxTagConst>(0,ncells), *this);
    Kokkos::fence();
    double time_const = timer.seconds() / ntrial / ncells;

    timer.reset();
    for (int i=0; i<ntrial; ++i)
      Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace,MomFluxTagConstTeam>(ncells,num_basis,32).set_scratch_size(0,Kokkos::PerThread(64*8*2)), *this);
    Kokkos::fence();
    double time_team = timer.seconds() / ntrial / ncells;

    printf("%5d   %9.3e   %9.3e   %9.3e   %9.3e   %9.3e\n",ncells,time_fad,time_fad_cont,time,time_const,time_team);

    if (do_check) {
      const double tol = 1e-14;
      check(residual_m_i, residual_m_i_fad, tol);
      check(residual_m_i, residual_m_i_fad_cont, tol);
      check(residual_m_i, residual_m_i_const, tol);
      check(residual_m_i, tresidual_m_i, tol);
    }
  }
};

template <typename ExecSpace>
void run(const int cell_begin, const int cell_end, const int cell_step,
         const int nbasis, const int npoint, const int ntrial, const bool check)
{
  const int fad_dim = 50;
  const int dim = 3;
  typedef DrekarTest<ExecSpace,dim,fad_dim> test_type;

  // The kernel allocates 2*N double's per warp on Cuda.  Approximate
  // the maximum number of warps as the maximum concurrency / 32.
  // Include a fudge factor of 1.2 since memory pool treats a block as full
  // once it reaches 80% capacity
  std::cout << "concurrency = " << ExecSpace::concurrency() << std::endl;
  const size_t block_size = fad_dim*sizeof(double);
  size_t nkernels = ExecSpace::concurrency()*2;
#if defined(KOKKOS_ENABLE_CUDA)
  if (std::is_same<ExecSpace, Kokkos::Cuda>::value)
    nkernels /= 32;
#endif
  size_t mem_pool_size =
    static_cast<size_t>(1.2*nkernels*block_size);
  const size_t superblock_size = std::max<size_t>(nkernels / 100, 1) * block_size;
    std::cout << "Memory pool size = " << mem_pool_size / (1024.0 * 1024.0)
              << " MB" << std::endl;
    ExecSpace exec_space;
    Sacado::createGlobalMemoryPool(exec_space, mem_pool_size,
        block_size, block_size, superblock_size);

#if defined(SACADO_VIEW_CUDA_HIERARCHICAL)
  std::cout << "hierarchical";
#elif defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
  std::cout << "hierarchical_dfad";
#else
  std::cout << "partitioned";
#endif
  std::cout << ", " << FadTypeName<typename test_type::FadType>::eval()
            << ":" << std::endl;

  printf("ncell      flat        hier      analytic     const        team\n");
  for(int i=cell_begin; i<=cell_end; i+=cell_step) {
    test_type test(i,nbasis,npoint);
    test.compute(ntrial, check);
  }

  Sacado::destroyGlobalMemoryPool(exec_space);
}

int main(int argc, char* argv[]) {
  bool success = true;
  try {

    // Set up command line options
    Teuchos::CommandLineProcessor clp(false);
    clp.setDocString("This program tests the speed of various forward mode AD implementations for simple Kokkos kernel");
#ifdef KOKKOS_ENABLE_SERIAL
    bool serial = 0;
    clp.setOption("serial", "no-serial", &serial, "Whether to run Serial");
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    int openmp = 0;
    clp.setOption("openmp", &openmp, "Number of OpenMP threads");
#endif
#ifdef KOKKOS_ENABLE_THREADS
    int threads = 0;
    clp.setOption("threads", &threads, "Number of pThreads threads");
#endif
#ifdef KOKKOS_ENABLE_CUDA
    bool cuda = 0;
    clp.setOption("cuda", "no-cuda", &cuda, "Whether to run CUDA");
#endif
    int numa = 0;
    clp.setOption("numa", &numa,
                  "Number of NUMA domains to use (set to 0 to use all NUMAs");
    int cores_per_numa = 0;
    clp.setOption("cores-per-numa", &cores_per_numa,
                  "Number of CPU cores per NUMA to use (set to 0 to use all cores)");
    bool print_config = false;
    clp.setOption("print-config", "no-print-config", &print_config,
                  "Whether to print Kokkos device configuration");
    int cell_begin = 100;
    clp.setOption("begin", &cell_begin, "Starting number of cells");
    int cell_end = 8000;
    clp.setOption("end", &cell_end, "Ending number of cells");
    int cell_step = 100;
    clp.setOption("step", &cell_step, "Cell increment");
    int nbasis = 8;
    clp.setOption("basis", &nbasis, "Number of basis functions");
    int npoint = 8;
    clp.setOption("point", &npoint, "Number of integration points");
    int ntrial = 5;
    clp.setOption("trial", &ntrial, "Number of trials");
    bool check = false;
    clp.setOption("check", "no-check", &check,
                  "Check correctness of results");

    // Parse options
    switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return 0;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return 1;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      break;
    }

    Kokkos::InitArguments init_args;
    init_args.num_threads = -1;
    #ifdef KOKKOS_ENABLE_OPENMP
      if(openmp) init_args.num_threads = openmp;
    #endif
    #ifdef KOKKOS_ENABLE_THREADS
      if(threads) init_args.num_threads = threads;
    #endif

    Kokkos::initialize(init_args);
    if (print_config)
      Kokkos::print_configuration(std::cout, true);

#ifdef KOKKOS_ENABLE_SERIAL
    if (serial) {
      using Kokkos::Serial;
      run<Serial>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_OPENMP
    if (openmp) {
      using Kokkos::OpenMP;
      run<OpenMP>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_THREADS
    if (threads) {
      using Kokkos::Threads;
      run<Threads>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif

#ifdef KOKKOS_ENABLE_CUDA
    if (cuda) {
      using Kokkos::Cuda;
      run<Cuda>(cell_begin, cell_end, cell_step, nbasis, npoint, ntrial, check);
    }
#endif
    Kokkos::finalize();
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return !success;
}
