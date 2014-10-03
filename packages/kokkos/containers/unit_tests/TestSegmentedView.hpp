/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_TEST_SEGMENTEDVIEW_HPP
#define KOKKOS_TEST_SEGMENTEDVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Core.hpp>
#include <Kokkos_SegmentedView.hpp>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

namespace Impl {

  template<class ViewType , class ExecutionSpace, int Rank = ViewType::Rank>
  struct GrowTest;

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 1> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+team_member.team_size());
      value += team_idx + team_member.team_rank();

      if((a.dimension_0()>team_idx+team_member.team_rank()) &&
         (a.dimension(0)>team_idx+team_member.team_rank()))
        a(team_idx+team_member.team_rank()) = team_idx+team_member.team_rank();
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 2> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        value += team_idx + team_member.team_rank() + 13*k;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++) {
          a(team_idx+ team_member.team_rank(),k) =
              team_idx+ team_member.team_rank() + 13*k;
        }
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 3> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          value += team_idx + team_member.team_rank() + 13*k + 3*l;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            a(team_idx+ team_member.team_rank(),k,l) =
                team_idx+ team_member.team_rank() + 13*k + 3*l;
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 4> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          for(int m=0;m<2;m++)
            value += team_idx + team_member.team_rank() + 13*k + 3*l + 7*m;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              a(team_idx+ team_member.team_rank(),k,l,m) =
                  team_idx+ team_member.team_rank() + 13*k + 3*l + 7*m;
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 5> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          for(int m=0;m<2;m++)
            for(int n=0;n<3;n++)
              value +=
                  team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                a(team_idx+ team_member.team_rank(),k,l,m,n) =
                  team_idx+ team_member.team_rank() + 13*k + 3*l + 7*m + 5*n;
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 6> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          for(int m=0;m<2;m++)
            for(int n=0;n<3;n++)
              for(int o=0;o<2;o++)
              value +=
                  team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o ;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                a(team_idx+ team_member.team_rank(),k,l,m,n,o) =
                    team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o ;
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 7> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      a.grow(team_member , team_idx+ team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          for(int m=0;m<2;m++)
            for(int n=0;n<3;n++)
              for(int o=0;o<2;o++)
                for(int p=0;p<4;p++)
              value +=
                  team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o + 15*p ;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                  for(int p=0;p<a.dimension_6();p++)
                a(team_idx+ team_member.team_rank(),k,l,m,n,o,p) =
                    team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o + 15*p ;
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct GrowTest<ViewType , ExecutionSpace , 8> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    GrowTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if(team_idx + team_member.team_size() > 100100 && team_idx + team_member.team_size() < 100110)
      printf("Size: %i %i %i %i | %i\n",team_member.league_rank(),team_member.team_rank(),team_member.league_size(),team_member.team_size(),team_idx + team_member.team_size());
      a.grow(team_member , team_idx + team_member.team_size());

      for(int k=0;k<7;k++)
        for(int l=0;l<3;l++)
          for(int m=0;m<2;m++)
            for(int n=0;n<3;n++)
              for(int o=0;o<2;o++)
                for(int p=0;p<4;p++)
                  for(int q=0;q<3;q++)
              value +=
                  team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o + 15*p + 17*q;

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                  for(int p=0;p<a.dimension_6();p++)
                    for(int q=0;q<a.dimension_7();q++)
                a(team_idx+ team_member.team_rank(),k,l,m,n,o,p,q) =
                    team_idx + team_member.team_rank() + 13*k + 3*l + 7*m + 5*n + 2*o + 15*p + 17*q;
      }
    }
  };

  template<class ViewType , class ExecutionSpace, int Rank = ViewType::Rank>
  struct VerifyTest;

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 1> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        value += a(team_idx+ team_member.team_rank());
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 2> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          value += a(team_idx+ team_member.team_rank(),k);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 3> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            value += a(team_idx+ team_member.team_rank(),k,l);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 4> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              value += a(team_idx+ team_member.team_rank(),k,l,m);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 5> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                value += a(team_idx+ team_member.team_rank(),k,l,m,n);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 6> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                  value += a(team_idx+ team_member.team_rank(),k,l,m,n,o);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 7> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                  for(int p=0;p<a.dimension_6();p++)
                    value += a(team_idx+ team_member.team_rank(),k,l,m,n,o,p);
      }
    }
  };

  template<class ViewType , class ExecutionSpace>
  struct VerifyTest<ViewType , ExecutionSpace , 8> {
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;
    typedef typename Policy::member_type team_type;
    typedef double value_type;

    ViewType a;

    VerifyTest(ViewType in):a(in) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (team_type team_member, double& value) const {
      unsigned int team_idx = team_member.league_rank() * team_member.team_size();

      if((a.dimension_0()>team_idx+ team_member.team_rank()) &&
         (a.dimension(0)>team_idx+ team_member.team_rank())) {
        for(int k=0;k<a.dimension_1();k++)
          for(int l=0;l<a.dimension_2();l++)
            for(int m=0;m<a.dimension_3();m++)
              for(int n=0;n<a.dimension_4();n++)
                for(int o=0;o<a.dimension_5();o++)
                  for(int p=0;p<a.dimension_6();p++)
                    for(int q=0;q<a.dimension_7();q++)
                      value += a(team_idx+ team_member.team_rank(),k,l,m,n,o,p,q);
      }
    }
  };

  template <typename Scalar, class ExecutionSpace>
  struct test_segmented_view
  {
    typedef test_segmented_view<Scalar,ExecutionSpace> self_type;

    typedef Scalar scalar_type;
    typedef ExecutionSpace execution_space;
    typedef Kokkos::TeamPolicy<execution_space> Policy;

    double result;
    double reference;

    template <class ViewType>
    void run_me(ViewType a, int max_length){
      reference = 0;
      result = 0;
      int team_size = execution_space::team_max();
#ifdef KOKKOS_HAVE_CUDA
      if(Kokkos::Impl::is_same<execution_space,Kokkos::Cuda>::value)
        team_size = 128;
      else
#endif
      if(team_size>4)
        team_size = 4;


      int nteams = max_length/team_size;

      Kokkos::parallel_reduce(Policy(nteams,team_size),GrowTest<ViewType,execution_space>(a),reference);
      Kokkos::parallel_reduce(Policy(nteams,team_size),VerifyTest<ViewType,execution_space>(a),result);
      Kokkos::fence();
    }


    test_segmented_view(unsigned int size,int rank)
    {
      reference = 0;
      result = 0;

      const int dim_1 = 7;
      const int dim_2 = 3;
      const int dim_3 = 2;
      const int dim_4 = 3;
      const int dim_5 = 2;
      const int dim_6 = 4;
      //const int dim_7 = 3;

      if(rank==1) {
        typedef Kokkos::SegmentedView<Scalar*,Kokkos::LayoutLeft,ExecutionSpace> rank1_view;
        run_me< rank1_view >(rank1_view("Rank1",128,size), size);
      }
      if(rank==2) {
        typedef Kokkos::SegmentedView<Scalar**,Kokkos::LayoutLeft,ExecutionSpace> rank2_view;
        run_me< rank2_view >(rank2_view("Rank2",128,size,dim_1), size);
      }
      if(rank==3) {
        typedef Kokkos::SegmentedView<Scalar*[7][3][2],Kokkos::LayoutRight,ExecutionSpace> rank3_view;
        run_me< rank3_view >(rank3_view("Rank3",128,size), size);
      }
      if(rank==4) {
        typedef Kokkos::SegmentedView<Scalar****,Kokkos::LayoutRight,ExecutionSpace> rank4_view;
        run_me< rank4_view >(rank4_view("Rank4",128,size,dim_1,dim_2,dim_3), size);
      }
      if(rank==5) {
        typedef Kokkos::SegmentedView<Scalar*[7][3][2][3],Kokkos::LayoutLeft,ExecutionSpace> rank5_view;
        run_me< rank5_view >(rank5_view("Rank5",128,size), size);
      }
      if(rank==6) {
        typedef Kokkos::SegmentedView<Scalar*****[2],Kokkos::LayoutRight,ExecutionSpace> rank6_view;
        run_me< rank6_view >(rank6_view("Rank6",128,size,dim_1,dim_2,dim_3,dim_4), size);
      }
      if(rank==7) {
        typedef Kokkos::SegmentedView<Scalar*******,Kokkos::LayoutLeft,ExecutionSpace> rank7_view;
        run_me< rank7_view >(rank7_view("Rank7",128,size,dim_1,dim_2,dim_3,dim_4,dim_5,dim_6), size);
      }
      if(rank==8) {
        typedef Kokkos::SegmentedView<Scalar*****[2][4][3],Kokkos::LayoutLeft,ExecutionSpace> rank8_view;
        run_me< rank8_view >(rank8_view("Rank8",128,size,dim_1,dim_2,dim_3,dim_4), size);
      }
    }

   };

} // namespace Impl




template <typename Scalar, class ExecutionSpace>
void test_segmented_view(unsigned int size)
{
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,1);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,2);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,3);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,4);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,5);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,6);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,7);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    Impl::test_segmented_view<Scalar,ExecutionSpace> test(size,8);
    ASSERT_EQ(test.reference,test.result);
  }
  {
    typedef Kokkos::SegmentedView<Scalar*****[2][4][3],Kokkos::LayoutLeft,ExecutionSpace> view_type;
    view_type a("A",128,size,7,3,2,1);
    double reference;
    int team_size = 1;
    int nteams = size/team_size;

    Kokkos::parallel_reduce(Kokkos::TeamPolicy<ExecutionSpace>(nteams,team_size),Impl::GrowTest<view_type,ExecutionSpace>(a),reference);
    size_t real_size = ((size+127)/128)*128;
    ASSERT_EQ(real_size,a.dimension_0());
    ASSERT_EQ(7,a.dimension_1());
    ASSERT_EQ(3,a.dimension_2());
    ASSERT_EQ(2,a.dimension_3());
    ASSERT_EQ(1,a.dimension_4());
    ASSERT_EQ(2,a.dimension_5());
    ASSERT_EQ(4,a.dimension_6());
    ASSERT_EQ(3,a.dimension_7());
    ASSERT_EQ(real_size,a.dimension(0));
    ASSERT_EQ(7,a.dimension(1));
    ASSERT_EQ(3,a.dimension(2));
    ASSERT_EQ(2,a.dimension(3));
    ASSERT_EQ(1,a.dimension(4));
    ASSERT_EQ(2,a.dimension(5));
    ASSERT_EQ(4,a.dimension(6));
    ASSERT_EQ(3,a.dimension(7));
    ASSERT_EQ(8,a.Rank);
  }
}


} // namespace Test

#endif //KOKKOS_TEST_UNORDERED_MAP_HPP
