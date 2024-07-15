// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//**********************************************************************
template<typename EvalT, typename Traits>
AllRanksEvaluator<EvalT,Traits>::
AllRanksEvaluator(const Teuchos::RCP<PHX::DataLayout>& dl1,
                  const Teuchos::RCP<PHX::DataLayout>& dl2,
                  const Teuchos::RCP<PHX::DataLayout>& dl3,
                  const Teuchos::RCP<PHX::DataLayout>& dl4,
                  const Teuchos::RCP<PHX::DataLayout>& dl5,
                  const Teuchos::RCP<PHX::DataLayout>& dl6) :
  f1("f1",dl1),
  f2("f2",dl2),
  f3("f3",dl3),
  f4("f4",dl4),
  f5("f5",dl5),
  f6("f6",dl6),
  x1("x1",dl1),
  x2("x2",dl2),
  x3("x3",dl3),
  x4("x4",dl4),
  x5("x5",dl5),
  x6("x6",dl6)
{
  this->addEvaluatedField(x1);
  this->addEvaluatedField(x2);
  this->addEvaluatedField(x3);
  this->addEvaluatedField(x4);
  this->addEvaluatedField(x5);
  this->addEvaluatedField(x6);
  this->addDependentField(f1);
  this->addDependentField(f2);
  this->addDependentField(f3);
  this->addDependentField(f4);
  this->addDependentField(f5);
  this->addDependentField(f6);
  this->setName("AllRanksEvaluator");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void AllRanksEvaluator<EvalT,Traits>::evaluateFields(typename Traits::EvalData workset_size)
{
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset_size,Kokkos::AUTO()),*this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void AllRanksEvaluator<EvalT,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int i = team.league_rank();

  Kokkos::single(Kokkos::PerTeam(team), [=] () {
      x1(i) = f1(i) * f1(i);
    });

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,static_cast<int>(x6.extent(1))), [&] (const int& j) {
      x2(i,j) = f2(i,j) * f2(i,j);
      for (int k = 0; k < static_cast<int>(x6.extent(2)); ++k) {
        x3(i,j,k) = f3(i,j,k) * f3(i,j,k);
        for (int l = 0; l < static_cast<int>(x6.extent(3)); ++l) {
          x4(i,j,k,l) = f4(i,j,k,l) * f4(i,j,k,l);
          for (int m = 0; m < static_cast<int>(x6.extent(4)); ++m) {
            x5(i,j,k,l,m) = f5(i,j,k,l,m) * f5(i,j,k,l,m);
            for (int n = 0; n < static_cast<int>(x6.extent(5)); ++n) {
              x6(i,j,k,l,m,n) = f6(i,j,k,l,m,n) * f6(i,j,k,l,m,n);
            }
          }
        }
      }
  });
}

//**********************************************************************
