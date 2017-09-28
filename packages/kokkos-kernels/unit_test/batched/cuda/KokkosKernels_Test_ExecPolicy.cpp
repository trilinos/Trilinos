/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <iomanip>

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

namespace KokkosKernels {

  namespace Test {

    // Assume that arrays are packed with VectorLength
#define VectorLength 8
#define ChunkPerThread VectorLength
#define ChunkPerTeam 96
#define N 100000

    enum { TEST_ADD = 0,
           TEST_MINUS = 1,
           TEST_MULT = 2,
           TEST_DIV = 3,
           TEST_UNARY_MINUS = 4 };

    struct TeamThreadTag {};
    struct ThreadVectorTag {};
    struct TeamThreadVectorTag {};
    struct ThreadVectorTeamTag {};
    
    template<typename ViewType, int TestID>
    struct Functor {
      ViewType _a, _b, _c;
      
      KOKKOS_INLINE_FUNCTION
      Functor(ViewType a, ViewType b, ViewType c)
        : _a(a), _b(b), _c(c) {}

      // range policy functor
      KOKKOS_INLINE_FUNCTION
      void operator()(const int i) const {
        switch (TestID) {
        case 0: _c(i) += (_a(i) + _b(i)); break;
        case 1: _c(i) += (_a(i) - _b(i)); break;
        case 2: _c(i) += (_a(i) * _b(i)); break;
        case 3: _c(i) += (_a(i) / _b(i)); break;
        case 4: _c(i) = -_c(i); break;
        }
      }
      
      // team policy functor: invoked with vector_length = 1
      template<typename TeamPolicyMemberType>
      KOKKOS_INLINE_FUNCTION
      void operator()(const TeamThreadTag&, const TeamPolicyMemberType &member) const {
        const int array_size = _c.dimension_0();
        const int ibegin = member.league_rank()*ChunkPerTeam;
        const int itmp = ibegin+ChunkPerTeam;
        const int iend = itmp < array_size ? itmp : array_size; 
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, ibegin, iend),
                             [&](const int i) {
			       operator()(i);
                             });
      }
      
      // vector policy functor: invoked with team_size = 1 
      template<typename TeamPolicyMemberType>
      KOKKOS_INLINE_FUNCTION
      void operator()(const ThreadVectorTag&, const TeamPolicyMemberType &member) const {
        const int array_size = _c.dimension_0();
        const int ibegin = member.league_rank()*ChunkPerTeam;
        const int itmp = ibegin+ChunkPerTeam;
        const int iend = itmp < array_size ? itmp : array_size; 
        
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, iend - ibegin),
                             [&](const int ii) {
                               const int i = ii + ibegin;
			       operator()(i);
                             });
      }

      // team-vector policy functor
      template<typename TeamPolicyMemberType>
      KOKKOS_INLINE_FUNCTION
      void operator()(const TeamThreadVectorTag&, const TeamPolicyMemberType &member) const {
        const int array_size = _c.dimension_0();

        const int ibegin = member.league_rank()*ChunkPerTeam;
        const int itmp = ibegin+ChunkPerTeam;
        const int iend = itmp < array_size ? itmp : array_size; 
        const int idiff = iend - ibegin;
        const int icnt = idiff/ChunkPerThread + (idiff%ChunkPerThread > 0);

        Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, icnt),
                             [&](const int idx_team) {
                               const int ibegin_thread = ibegin + idx_team*ChunkPerThread;
                               const int itmp_thread = ibegin_thread + ChunkPerThread;
                               const int iend_thread = itmp_thread < iend ? itmp_thread : iend; 
                               const int icnt_thread = iend_thread - ibegin_thread;

                               Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, icnt_thread),
                                                    [&](const int idx_thread) {
                                                      const int i = idx_thread + ibegin_thread;
						      operator()(i);
                                                    });
                             });
      }

      // vector-team policy functor
      // 
      template<typename TeamPolicyMemberType>
      KOKKOS_INLINE_FUNCTION
      void operator()(const ThreadVectorTeamTag&, const TeamPolicyMemberType &member) const {
        const int array_size = _c.dimension_0();

        const int ibegin = member.league_rank()*ChunkPerTeam;
        const int itmp = ibegin+ChunkPerTeam;
        const int iend = itmp < array_size ? itmp : array_size; 

        const int idiff = iend - ibegin;
        const int icnt = idiff/ChunkPerThread + (idiff%ChunkPerThread > 0);

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, ChunkPerThread),
                             [&](const int idx_thread) {
                               
                               Kokkos::parallel_for(Kokkos::TeamThreadRange(member, 0, icnt),
                                                    [&](const int ii) {
                                                      const int i = ibegin + ii*ChunkPerThread + idx_thread;
						      operator()(i);
                                                    });
                             });
      }

    };
    
    template<typename DeviceSpaceType, typename ValueType, int TestID>
    void ExecPolicy() {
      typedef Kokkos::DefaultHostExecutionSpace HostSpaceType;

      const int iter_begin = -10, iter_end = 100;
      Kokkos::Impl::Timer timer;

      {
        typedef Kokkos::View<ValueType*,HostSpaceType> HostViewType;
        HostViewType 
          a_host("a_host", N*VectorLength), 
          b_host("b_host", N*VectorLength), 
          c_host("c_host", N*VectorLength),
          cref_host("cref_host", N*VectorLength);

        for (int k=0;k<N*VectorLength;++k) {
          const int 
            i = k/VectorLength,
            j = k%VectorLength;
          a_host(k) = j + 1;
          b_host(k) = i + 1;
          c_host(k) = i*j;
        }
        
        auto a  = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), a_host);
        auto b  = Kokkos::create_mirror_view(typename DeviceSpaceType::memory_space(), b_host);

        typedef Kokkos::View<ValueType*,DeviceSpaceType> DeviceViewType;
        DeviceViewType c("c", N*VectorLength);
        
        Kokkos::deep_copy(a,  a_host);
        Kokkos::deep_copy(b,  b_host);
        Kokkos::deep_copy(c,  c_host);

        ///
        /// RangePolicy
        ///
        {
          double t = 0;
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::RangePolicy<DeviceSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, N*VectorLength);
            Kokkos::parallel_for( policy, Functor<DeviceViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          Kokkos::deep_copy(cref_host, c);
          std::cout << "test = " << std::setw(3) << TestID 
                    << " time(range policy)      = " << std::scientific << (t/iter_end) 
                    << std::endl;
        }

        ///
        /// TeamPolicy (vector size = 1)
        ///
        {
          Kokkos::deep_copy(c, c_host);
          
          double t = 0;
          const int lsize = c.dimension_0()/ChunkPerTeam + (c.dimension_0()%ChunkPerTeam > 0);
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::TeamPolicy<DeviceSpaceType,TeamThreadTag,Kokkos::Schedule<Kokkos::Static> > policy(lsize, Kokkos::AUTO, 1);
            Kokkos::parallel_for( policy, Functor<DeviceViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          auto csol_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol_host, c);

          double sum = 0;
          for (int i=0;i<cref_host.dimension_0();++i) 
            sum += std::abs(cref_host(i) - csol_host(i));
          std::cout << "test = " << std::setw(3) << TestID 
                    << " time(team policy)       = " << std::scientific << (t/iter_end) 
                    << " diff = " << sum << std::endl;
        }

        ///
        /// TeamPolicy (team size = 1)
        ///
        {
          Kokkos::deep_copy(c, c_host);
          
          double t = 0;
          const int lsize = c.dimension_0()/ChunkPerTeam + (c.dimension_0()%ChunkPerTeam > 0);
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::TeamPolicy<DeviceSpaceType,ThreadVectorTag,Kokkos::Schedule<Kokkos::Static> > policy(lsize, 1, VectorLength);
            Kokkos::parallel_for( policy, Functor<DeviceViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          auto csol_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol_host, c);

          double sum = 0;
          for (int i=0;i<cref_host.dimension_0();++i)
            sum += std::abs(cref_host(i) - csol_host(i));
          std::cout << "test = " << std::setw(3) << TestID 
                    << " time(vector policy)     = " << std::scientific << (t/iter_end) 
                    << " diff = " << sum << std::endl;
        }

        ///
        /// TeamPolicy (team - vector)
        ///
        {
          Kokkos::deep_copy(c, c_host);
          
          double t = 0;
          const int lsize = c.dimension_0()/ChunkPerTeam + (c.dimension_0()%ChunkPerTeam > 0);
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::TeamPolicy<DeviceSpaceType,TeamThreadVectorTag,Kokkos::Schedule<Kokkos::Static> > policy(lsize, Kokkos::AUTO, VectorLength);
            Kokkos::parallel_for( policy, Functor<DeviceViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          auto csol_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol_host, c);

          double sum = 0;
          for (int i=0;i<cref_host.dimension_0();++i)
            sum += std::abs(cref_host(i) - csol_host(i));
          std::cout << "test = " << std::setw(3) << TestID 
                    << " time(team-vector policy)= " << std::scientific << (t/iter_end) 
                    << " diff = " << sum << std::endl;
        }
        
        ///
        /// TeamPolicy (vector - team)
        ///
        {
          Kokkos::deep_copy(c, c_host);
          
          double t = 0;
          const int lsize = c.dimension_0()/ChunkPerTeam + (c.dimension_0()%ChunkPerTeam > 0);
          for (int iter=iter_begin;iter<iter_end;++iter) {
            DeviceSpaceType::fence();
            timer.reset();
            
            Kokkos::TeamPolicy<DeviceSpaceType,ThreadVectorTeamTag,Kokkos::Schedule<Kokkos::Static> > policy(lsize, Kokkos::AUTO, VectorLength);
            Kokkos::parallel_for( policy, Functor<DeviceViewType,TestID>(a, b, c) );
            
            DeviceSpaceType::fence();
            t += (iter >= 0)*timer.seconds();
          }
          auto csol_host = Kokkos::create_mirror_view(typename HostSpaceType::memory_space(), c);
          Kokkos::deep_copy(csol_host, c);

          double sum = 0;
          for (int i=0;i<cref_host.dimension_0();++i)
            sum += std::abs(cref_host(i) - csol_host(i));
          std::cout << "test = " << std::setw(3) << TestID 
                    << " time(vector-team policy)= " << std::scientific << (t/iter_end) 
                    << " diff = " << sum << std::endl;
        }
        
      }

    }
  }
}

int main(int argc, char *argv[]) {

  Kokkos::initialize();

  std::cout << "  VectorLength = " <<  VectorLength 
	    << "  ChunkPerThread = " << ChunkPerThread
	    << "  ChunkPerTeam = " << ChunkPerTeam
	    << "  N = " << N
	    << "\n";

  std::cout << "\n\n Testing Serial Range, Team, Vector Policy double \n";
#if defined(KOKKOS_HAVE_SERIAL)
  KokkosKernels::Test::ExecPolicy<Kokkos::Serial,double,0>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Serial,double,1>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Serial,double,2>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Serial,double,3>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Serial,double,4>();
#else
  std::cout << "Kokkos::Serial is not enabled\n";
#endif

  std::cout << "\n\n Testing OpenMP Range, Team, Vector Policy double \n";
#if defined(KOKKOS_HAVE_OPENMP)
  KokkosKernels::Test::ExecPolicy<Kokkos::OpenMP,double,0>();
  KokkosKernels::Test::ExecPolicy<Kokkos::OpenMP,double,1>();
  KokkosKernels::Test::ExecPolicy<Kokkos::OpenMP,double,2>();
  KokkosKernels::Test::ExecPolicy<Kokkos::OpenMP,double,3>();
  KokkosKernels::Test::ExecPolicy<Kokkos::OpenMP,double,4>();
#else
  std::cout << "Kokkos::OpenMP is not enabled\n";
#endif

  std::cout << "\n\n Testing Cuda Range, Team, Vector Policy double \n";
#if defined(KOKKOS_HAVE_CUDA)
  KokkosKernels::Test::ExecPolicy<Kokkos::Cuda,double,0>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Cuda,double,1>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Cuda,double,2>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Cuda,double,3>();
  KokkosKernels::Test::ExecPolicy<Kokkos::Cuda,double,4>();
#else 
  std::cout << "Kokkos::Cuda is not enabled\n";
#endif
  Kokkos::finalize();

  return 0;
}

