#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <limits>

#ifdef KOKKOS_HAVE_CUDA
  typedef Kokkos::Cuda Device;
  typedef Kokkos::CudaSpace MemSpace;
#else
  typedef Kokkos::Serial Device;
  typedef Kokkos::HostSpace MemSpace;
#endif

template <typename T>
struct MinOp
{
    KOKKOS_FUNCTION
    void operator()(volatile T& update, volatile const T& input) const
    {
        update = update < input ? update : input;
    }
};

typedef double value_type;
typedef Kokkos::View<value_type*, Kokkos::LayoutRight, MemSpace> view_type;

template<typename ReductionOp>
struct ViewAccessFunctor
{
  KOKKOS_FUNCTION
  ViewAccessFunctor(const view_type& view, int team_offset) : v(view), offset(team_offset) {}

  KOKKOS_FUNCTION
  void operator()(int i, value_type& update) const
  {
    ReductionOp()(update,v(offset + i));
  }

private:
  const view_type& v;
  int offset;
};

template <typename ReductionOp>
struct ReductionTeamFunctor
{
    struct JoinOp {
    public:
        typedef JoinOp reducer;
        typedef double value_type;

    private:
        value_type* value;

    public:
        KOKKOS_FUNCTION
        JoinOp (value_type& value_) : value (&value_) {}

        KOKKOS_FUNCTION
        void join (value_type& dest, const value_type& src) const {
            ReductionOp() (dest, src);
        }

        KOKKOS_FUNCTION
        void join (volatile value_type& dest, const volatile value_type& src) const {
            ReductionOp() (dest, src);
        }

        KOKKOS_FUNCTION
        void init (value_type& val) const {
            ReductionOp()(val, *value);
        }

        KOKKOS_FUNCTION
        value_type& reference() const {
          return *value;
        }

        KOKKOS_FUNCTION
        view_type view () const {
            return view_type (value);
        }
    };

    KOKKOS_FUNCTION
    ReductionTeamFunctor(const view_type v, int len, value_type initVal) : view(v), inner_length(len), initialValue(initVal) { }

    KOKKOS_FUNCTION
    void init(value_type &update) const
    {
        update = initialValue;
    }

    KOKKOS_FUNCTION
    void join(volatile value_type& update, volatile const value_type& input) const
    {
        ReductionOp()(update, input);
    }

    typedef typename Kokkos::TeamPolicy<Device, Kokkos::Schedule<Kokkos::Dynamic> >::member_type TeamHandleType;

    KOKKOS_FUNCTION
    void operator()(const TeamHandleType& team, value_type& update) const
    {
        const int teamIndex = team.league_rank();
        unsigned numElements = inner_length;
        value_type localUpdate = initialValue;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, 0u, numElements),
                                ViewAccessFunctor<MinOp<double> >(view, teamIndex*inner_length),
                                JoinOp(localUpdate));
        Kokkos::single(Kokkos::PerTeam(team), [&](){join(update, localUpdate);});
    }

private:
    view_type view;
    int inner_length;
    value_type initialValue;
};

double get_min_value()
{
    const int num_teams = 5;
    const int inner_length = 5;
    const int total_length = num_teams * inner_length;

    Kokkos::View<double*,Kokkos::LayoutRight,MemSpace> data("data", num_teams*inner_length);

    Kokkos::parallel_for(total_length, KOKKOS_LAMBDA(const int& i) {
        data(i) = i+1;
    });

    double reduction_result = std::numeric_limits<double>::max();
    ReductionTeamFunctor<MinOp<double> > teamFunctor(data, inner_length, reduction_result);

    Kokkos::parallel_reduce(Kokkos::TeamPolicy<Device>(num_teams, Kokkos::AUTO), teamFunctor, reduction_result);

    //min value should be 1 !!
    return reduction_result;
}

TEST(KokkosReduce, get_min_value)
{
    EXPECT_EQ(1, get_min_value());
}

