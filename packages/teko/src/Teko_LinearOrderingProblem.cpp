// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_HeuristicBlockPermutation.hpp"
#include <limits>
#include <vector>
#include <algorithm>

// Linear ordering problem solver adapted from:
//     CHARON, I. AND HUDRY, O.: A branch-and-bound algorithm to solve the linear ordering
//     problem for weighted tournaments, Discrete Applied Mathematics 156 (2006), 2097-2116
namespace Teko::LinearOrdering {
using ViewType = BlockNormsViewType;

namespace {

ViewType convert_to_normal_form(const ViewType &blockNorms) {
  ViewType normalForm("normalForm", blockNorms.extent(0), blockNorms.extent(1));
  for (auto i = 0U; i < blockNorms.extent(0); ++i) {
    for (auto j = 0U; j < blockNorms.extent(1); ++j) {
      const auto v     = blockNorms(i, j);
      normalForm(i, j) = v * v;
      if (i == j) normalForm(i, i) = 0.0;
    }
  }

  for (auto i = 0U; i < blockNorms.extent(0); ++i) {
    for (auto j = 0U; j <= i; ++j) {
      const auto hij   = normalForm(i, j);
      const auto hji   = normalForm(j, i);
      normalForm(i, j) = hij - hji;
      normalForm(j, i) = -normalForm(i, j);
    }
  }

  return normalForm;
}

struct Circuit {
  unsigned int vertex1, vertex2, vertex3;
  double value;
};

struct Node {
  int vertex;
  double value;
  double LowerBound;
  std::shared_ptr<Node> child;
  std::shared_ptr<Node> brother;
};

std::shared_ptr<Node> create(const std::shared_ptr<Node> &current, bool child, int vertex) {
  auto pointer        = std::make_shared<Node>();
  pointer->vertex     = vertex;
  pointer->value      = -1;
  pointer->LowerBound = -1;
  pointer->child      = nullptr;
  if (child) {
    pointer->brother = current->child;
    current->child   = pointer;
  } else {
    pointer->brother = current->brother;
    current->brother = pointer;
  }
  return pointer;
}

auto init_head_node() {
  auto head     = std::make_shared<Node>();
  head->vertex  = -1;
  head->value   = -1;
  head->child   = nullptr;
  head->brother = nullptr;
  return head;
}

void process_current_node(std::shared_ptr<Node> &current, int j, bool &IsInTree) {
  if (current->child == nullptr) {
    current  = create(current, false, j);
    IsInTree = false;
    return;
  }

  if (current->child->vertex > j) {
    current  = create(current, false, j);
    IsInTree = false;
    return;
  }

  current = current->child;
  if (current->vertex == j) return;

  while (current->brother != nullptr) {
    if (current->brother->vertex <= j) {
      current = current->brother;
      if (current->vertex == j) return;
    } else {
      break;
    }
  }
  if (current->vertex != j) {
    current  = create(current, true, j);
    IsInTree = false;
  }
}

void change(int oldPlace, int newPlace, int *order, int *place) {
  if (newPlace == oldPlace) return;

  const auto movingVertex = order[oldPlace];
  if (newPlace < oldPlace) {
    for (int i = oldPlace; i > newPlace; i--) {
      order[i]        = order[i - 1];
      place[order[i]] = i;
    }
  } else {
    for (int i = oldPlace; i < newPlace; i++) {
      order[i]        = order[i + 1];
      place[order[i]] = i;
    }
  }

  order[newPlace]     = movingVertex;
  place[movingVertex] = newPlace;
}

double adjust_ro_based_on_lagrangian(double L, double oldL, double ro) {
  if ((L - oldL) / L > 0.1) {
    ro *= 1.1;
    return ro;
  }

  if ((L - oldL) / L < -0.95) {
    ro *= 0.8;
    return ro;
  }

  ro *= 1 + (L - oldL) / L;
  return ro;
}

enum class TerminationReason { NONE, OPTIMAL_SOLUTION, MAX_WALL_TIME };

std::string to_string(TerminationReason reason) {
  switch (reason) {
    case TerminationReason::NONE: return "NONE";
    case TerminationReason::OPTIMAL_SOLUTION: return "finding optimal solution";
    case TerminationReason::MAX_WALL_TIME: return "hitting max walltime";
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidArgument,
                                 "Unknown termination reason.");
      return "";
  }
}

struct LinearOrderingSolver {
 public:
  LinearOrderingSolver(const ViewType &m_tournament, double m_tMaxWalltime,
                       const Teuchos::RCP<const Teuchos::Comm<int>> &m_communicator)
      : N(static_cast<int>(m_tournament.extent(0))),
        communicator(m_communicator),
        tMaxWalltime(m_tMaxWalltime),
        tournament(m_tournament),
        diffLambda("diffLambda", N, N, N),
        saveCoeff("saveCoeff", N, N),
        benefit("benefit", N, N),
        X("X", N, N),
        currentOrder(N),
        bestOrder(N),
        score(N),
        place(N),
        isOutBegSec(N),
        sortedVertices(N),
        outBegSecBenefit(N) {
    std::iota(currentOrder.begin(), currentOrder.end(), 0);
    std::iota(bestOrder.begin(), bestOrder.end(), 0);
    std::fill(isOutBegSec.begin(), isOutBegSec.end(), true);
    bound = order_value(bestOrder.data(), N);

    if (N % 2) {
      MaxNbCircuits = (N * N * N - N) / 24;
    } else {
      MaxNbCircuits = (N * N * N - 4 * N) / 24;
    }
    circuits.resize(MaxNbCircuits);

    init_branch_and_bound();
  }

  void solve() {
    // currently, only rank 0 performs the search as this process is difficult to parallelize
    if (communicator->getRank() == 0) {
      terminationReason = TerminationReason::OPTIMAL_SOLUTION;
      branch_and_bound(0, 0);
    }

    int reason         = static_cast<int>(terminationReason);
    const int rootRank = 0;
    Teuchos::broadcast(*communicator, rootRank, N, bestOrder.data());
    Teuchos::broadcast(*communicator, rootRank, 1, &reason);
    terminationReason = static_cast<TerminationReason>(reason);
  }

  [[nodiscard]] const auto &order() const { return bestOrder; }
  [[nodiscard]] auto termination_reason() const { return terminationReason; }

 private:
  static constexpr int NbIterMax1 = 1000;
  static constexpr int NbIterMax2 = 10;
  TerminationReason terminationReason{TerminationReason::NONE};

  Kokkos::Timer timer{};

  int N;
  const Teuchos::RCP<const Teuchos::Comm<int>> &communicator;
  double tMaxWalltime;
  ViewType tournament;
  Kokkos::View<double ***, Kokkos::HostSpace> diffLambda;
  Kokkos::View<double **, Kokkos::HostSpace> saveCoeff;
  Kokkos::View<double **, Kokkos::HostSpace> benefit{};
  Kokkos::View<int **, Kokkos::HostSpace> X;
  std::vector<int> currentOrder;
  std::vector<int> bestOrder;
  std::vector<int> score;
  std::vector<int> place;
  double bound = std::numeric_limits<double>::max();
  std::vector<bool> isOutBegSec;
  bool firstRelaxation{true};
  std::vector<int> sortedVertices;
  double positiveTriangle{};
  double ro1{}, ro2{}, ro{};
  int NbCircuits{};
  std::vector<Circuit> circuits;
  int MaxNbCircuits{};
  std::vector<double> outBegSecBenefit;
  std::vector<std::vector<int>> toExam;
  std::shared_ptr<Node> currentNode{};
  std::shared_ptr<Node> head = nullptr;

  std::pair<double, int> best_place(int oldPlace, int localN, int *order) {
    double minimum = 0;
    double delta   = 0;
    int vertex     = order[oldPlace];

    int newPlaceAddr = oldPlace;
    for (int i = oldPlace - 1; i >= 0; i--) {
      delta += tournament(order[i], vertex);
      if (delta < minimum) {
        minimum      = delta;
        newPlaceAddr = i;
      }
    }
    delta = 0;
    for (int i = oldPlace + 1; i < localN; i++) {
      delta += tournament(vertex, order[i]);
      if (delta < minimum) {
        minimum      = delta;
        newPlaceAddr = i;
      }
    }
    return std::make_pair(minimum, newPlaceAddr);
  }

  double descent(int localN, int *order) {
    for (int i = 0; i < localN; i++) {
      place[order[i]] = i;
    }

    int index       = 0;
    int index1      = 0;
    int NbUnchanged = 0;
    do {
      const auto [delta, index2] = best_place(index1, localN, order);
      if (delta < 0) {
        change(index1, index2, order, place.data());
        NbUnchanged = 0;
      } else
        NbUnchanged++;
      index             = (index + 1) % localN;
      const auto vertex = order[index];
      index1            = place[vertex];
    } while (NbUnchanged < localN);

    return order_value(order, localN);
  }

  bool relax(int *LocalVertices, int localN, double LocalBound, double LocalbestOrder_value) {
    double L                 = 0;
    double gain              = 0.0;
    double bestL             = 0;
    double currentOrderValue = 0.0;
    bool prepareGain         = true;
    double OldML             = LocalBound;
    double ML                = 0.0;
    double rohigh            = ro1;
    double min               = 0.05;

    enum class TerminationStatus { Converged, Unconverged, MaxIterations };
    TerminationStatus end = TerminationStatus::Unconverged;

    NbCircuits    = 0;
    int NbIterMax = firstRelaxation ? NbIterMax1 : NbIterMax2;
    ro            = firstRelaxation ? ro1 : ro2;

    init_one_relax(LocalVertices, localN);
    double oldL = 0.0;

    ML = L = dual_function(localN);
    if (bestL >= LocalBound) end = TerminationStatus::Converged;

    int iterNum = 0;
    while (end == TerminationStatus::Unconverged) {
      iterNum++;
      if (update_diff_lambda(localN, L, LocalbestOrder_value, prepareGain)) {
        currentOrderValue = positiveTriangle;
        for (int i = 0; i < localN - 1; i++) {
          for (int j = i + 1; j < localN; j++) {
            if (X(i, j)) currentOrderValue += tournament(LocalVertices[j], LocalVertices[i]);
          }
        }

        if (currentOrderValue < LocalBound) {
          bound += currentOrderValue - LocalBound;
          LocalBound           = currentOrderValue;
          LocalbestOrder_value = currentOrderValue;
          save_order(LocalVertices, localN);
          if (L >= LocalBound) end = TerminationStatus::Converged;
        }
      } else if (prepareGain) {
        gain = compute_gain(localN);
        L += gain;
        if (L > bestL) bestL = L;
        if (bestL >= LocalBound) {
          end = TerminationStatus::Converged;
          break;
        }
        prepareGain = false;
        gain        = 0;
      }

      L     = dual_function(localN);
      ML    = std::max(ML, L);
      bestL = std::max(bestL, L);
      if (bestL >= LocalBound) end = TerminationStatus::Converged;

      if (firstRelaxation) {
        if (L < 0) {
          ro *= 0.5;
          ro2 *= 0.99;
          Kokkos::deep_copy(diffLambda, 0.0);
          oldL = 0;
        }
        if (iterNum > NbIterMax - 500) ro *= 0.999;
        ro = adjust_ro_based_on_lagrangian(L, oldL, ro);
        if (iterNum == 400) {
          OldML = bestL;
        } else if ((iterNum >= 500) && (!(iterNum % 100))) {
          if ((ML - OldML) / (LocalBound - OldML) > min) {
            NbIterMax += 100;
          }
          OldML = ML;
          if (ro > rohigh) rohigh = ro;
        }
      } else {
        if (L < 0) {
          ro *= 0.5;
          ro2 *= 0.99;
        } else if (L < oldL * 0.98) {
          ro *= 0.9;
          ro2 *= 0.99;
        }
        ro *= 0.95;
      }

      oldL = L;
      if (iterNum >= NbIterMax) end = TerminationStatus::MaxIterations;
    }

    if (firstRelaxation) {
      ro2 = rohigh * 0.8;
    }
    firstRelaxation = false;
    set_lower_bound(bestL);

    return end == TerminationStatus::Converged;
  }

  double dual_function(int localN) {
    double L = positiveTriangle;
    for (int i = 0; i < localN - 1; i++) {
      auto Xi        = Kokkos::subview(X, i, Kokkos::ALL);
      const auto SVi = sortedVertices[i];
      for (int j = i + 1; j < localN; j++) {
        const auto SVj = sortedVertices[j];
        auto coeff     = tournament(SVj, SVi);
        for (int k = j + 1; k < localN; k++) {
          auto diff = diffLambda(SVi, SVj, sortedVertices[k]);
          if (diff > 0) L -= diff;
          coeff += diff;
        }
        for (int k = 0; k < i; k++) coeff += diffLambda(sortedVertices[k], SVi, SVj);
        for (int k = i + 1; k < j; k++) coeff -= diffLambda(SVi, sortedVertices[k], SVj);
        if (coeff < 0) {
          L += coeff;
          Xi[j] = 1;
        } else
          Xi[j] = 0;
        saveCoeff(i, j) = coeff;
      }
    }
    return L;
  }

  bool update_diff_lambda(int localN, double L, double LocalbestOrder_value, bool PrepareGain) {
    double coij     = 0;
    int S           = 0;
    double D        = 0;
    bool transitive = true;
    double gain     = 0;

    const auto simpleStep =
        (ro * (LocalbestOrder_value - L)) / (localN * (localN - 1) * (localN - 2));
    const auto doubleStep    = 2 * simpleStep;
    const auto tripleStep    = 3 * simpleStep;
    const auto negSimpleStep = -simpleStep;
    const auto negDoubleStep = -doubleStep;
    const auto negTripleStep = -tripleStep;
    NbCircuits               = 0;

    auto updateDiffij = [&](double &Di, double &diffijOk, int Si) {
      switch (Si) {
        case 2:
          if (Di > 0) {
            diffijOk += simpleStep;
          } else if (Di > negDoubleStep) {
            diffijOk = simpleStep;
          } else {
            diffijOk += tripleStep;
          }
          break;
        case 1:
          if (Di < negSimpleStep) {
            diffijOk += simpleStep;
          } else if (Di < 0) {
            diffijOk = 0;
          }
          break;
        case 0:
          if (Di > simpleStep) {
            diffijOk += negSimpleStep;
          } else if (Di > 0) {
            diffijOk = 0;
          }
          break;
        case -1:
          if (Di > doubleStep) {
            diffijOk += negTripleStep;
          } else if (Di > 0) {
            diffijOk = negSimpleStep;
          } else {
            diffijOk += negSimpleStep;
          }
          break;
      }
    };

    auto handleCircuit = [&](int i, int j, int k, int Si, double coi_j, double coj_k,
                             double coi_k) {
      transitive = false;

      Circuit newCircuit{};
      newCircuit.vertex1 = i;
      newCircuit.vertex2 = j;
      newCircuit.vertex3 = k;
      if (Si == -1) {
        gain = std::min({coi_j, coj_k, -coi_k});
      } else if (Si == 2) {
        gain = std::min({-coi_j, -coj_k, coi_k});
      }
      newCircuit.value       = gain;
      circuits[NbCircuits++] = newCircuit;
    };

    for (int i = 0; i < localN - 2; i++) {
      auto diffi = Kokkos::subview(diffLambda, sortedVertices[i], Kokkos::ALL, Kokkos::ALL);
      auto Xi    = Kokkos::subview(X, i, Kokkos::ALL);
      auto coi   = Kokkos::subview(saveCoeff, i, Kokkos::ALL);
      for (int j = i + 1; j < localN - 1; j++) {
        auto Xj     = Kokkos::subview(X, j, Kokkos::ALL);
        auto coj    = Kokkos::subview(saveCoeff, j, Kokkos::ALL);
        auto diffij = Kokkos::subview(diffi, sortedVertices[j], Kokkos::ALL);
        coij        = PrepareGain ? coi[j] : 0;
        for (int k = j + 1; k < localN; k++) {
          const auto Ok = sortedVertices[k];
          S             = Xi[j] + Xj[k] - Xi[k];
          if (PrepareGain) {
            if (S == -1 || S == 2) {
              handleCircuit(i, j, k, S, coij, coj[k], coi[k]);
            }
          } else {
            if (S == -1 || S == 2) {
              transitive = false;
            }
          }
          D = diffij[Ok];
          updateDiffij(D, diffij[Ok], S);
        }
      }
    }
    return transitive;
  }

  void save_order(int *localVertices, int localN) {
    std::fill_n(score.begin(), localN, 0);
    for (int i = 0; i < localN - 1; i++) {
      for (int j = i + 1; j < localN; j++) {
        if (X(i, j)) {
          score[i]++;
        } else {
          score[j]++;
        }
      }
    }

    for (int i = 0; i < localN; i++) {
      int index = 0;
      while (score[index] != localN - i - 1) index++;
      localVertices[i] = sortedVertices[index];
    }

    std::copy_n(localVertices - N + localN, N, bestOrder.data());
  }

  double compute_gain(int localN) {
    std::sort(circuits.begin(), circuits.begin() + NbCircuits,
              [](const Circuit &a, const Circuit &b) { return a.value < b.value; });
    double gain = 0;

    const int sentinel = 10;  // entries in X take on values of 0, 1, 10, or 11 (with the latter two
                              // being used to mark whether to compute the gain)

    for (int i = 0; i < NbCircuits; i++) {
      auto ACircuit = circuits[i];
      if ((X(ACircuit.vertex1, ACircuit.vertex2) < sentinel) &&
          (X(ACircuit.vertex2, ACircuit.vertex3) < sentinel) &&
          (X(ACircuit.vertex1, ACircuit.vertex3) < sentinel)) {
        gain += ACircuit.value;
        X(ACircuit.vertex1, ACircuit.vertex2) += sentinel;
        X(ACircuit.vertex2, ACircuit.vertex3) += sentinel;
        X(ACircuit.vertex1, ACircuit.vertex3) += sentinel;
      }
    }

    for (int i = 0; i < localN - 1; i++) {
      for (int j = i + 1; j < localN; j++) {
        if (X(i, j) >= sentinel) X(i, j) -= sentinel;
      }
    }

    return gain;
  }

  void init_one_relax(int * /*LocalVertices*/, int localN) {
    positiveTriangle = 0;
    int k            = 0;

    for (int i = 0; i < N; i++) {
      if (isOutBegSec[i]) sortedVertices[k++] = i;
    }

    for (int i = 0; i < localN - 1; i++) {
      for (int j = i + 1; j < localN; j++) {
        if (tournament(sortedVertices[i], sortedVertices[j]) > 0) {
          positiveTriangle += tournament(sortedVertices[i], sortedVertices[j]);
        }
      }
    }
  }

  double order_value(int *order, int q) const {
    double weight = 0;
    for (int i = 0; i < q - 1; i++) {
      for (int j = i + 1; j < q; j++) {
        if (tournament(order[i], order[j]) < 0) {
          weight += tournament(order[j], order[i]);
        }
      }
    }
    return weight;
  }

  bool is_valid(int vertex, int index) const {
    /*test named "ham"*/
    if ((index >= 1) && (tournament(currentOrder[index - 1], vertex) < 0)) return false;

    /* special case of OSmoves, when we consider the move of "vertex" at
       the end of the total order*/
    if (outBegSecBenefit[vertex] < 0) return false;

    double CurrentBenefit = 0;
    for (int j = index - 1; j >= 0; j--) {
      /* special case of Smoves:  we consider the move of "vertex"
         before the vertex which is now at the index j*/
      CurrentBenefit += tournament(vertex, currentOrder[j]);
      if (CurrentBenefit > 0) {
        return false;
      } else if (CurrentBenefit == 0)
        if (currentOrder[j] > vertex) return false;
    }

    for (int j = 0; j < index - 1; j++) {
      /* special case of OSmoves: we consider the move of the vertex
         which is now at the index j after the beginning section obtained
         by adding "vertex" after the current beginning section*/
      CurrentBenefit = benefit(j, index - 1) + tournament(currentOrder[j], vertex);
      if (CurrentBenefit < 0) {
        return false;
      } else if (CurrentBenefit == 0)
        if (currentOrder[j] > currentOrder[j + 1]) return false;
    }

    /* special case of lex: if the weight between the last vertex
       of the current beginning section and "vertex" is equal to 0,
       and if the "name" of the last vertex of the current beginning
       section is greater (for the lexicographic order) than the "name"
       of "vertex", we do not study the beginning section obtained  by
       adding "vertex" after the current beginning section */
    if ((index >= 1) && (tournament(currentOrder[index - 1], vertex) == 0) &&
        (vertex < currentOrder[index - 1]))
      return false;

    /*test OSmoves, only performed for intervals having no more than
     3 vertices; for more vertices, we observed that it nearly never cut*/
    CurrentBenefit = outBegSecBenefit[vertex];
    for (int i = index - 1; (i >= 0) && (i >= index - 3); i--) {
      CurrentBenefit += outBegSecBenefit[currentOrder[i]] - tournament(currentOrder[i], vertex);
      if (CurrentBenefit < 0) return false;
    }

    /*test Smoves, only performed for intervals having no more than
     3 vertices; for more vertices, we observed that it nearly never cut*/
    for (int j = 1; j < index - 1; j++) {
      CurrentBenefit = benefit(j, index - 1) + tournament(currentOrder[j], vertex);
      for (int i = j - 1; i >= 0; i--) {
        if ((j - i > 3) && (index - j > 4)) break;
        CurrentBenefit +=
            benefit(i, index - 1) + tournament(currentOrder[i], vertex) - benefit(i, j);
        if (CurrentBenefit < 0) {
          return false;
        } else if (CurrentBenefit == 0) {
          /* lexicographic test */
          if (currentOrder[i] > currentOrder[j + 1]) return false;
        }
      }
    }
    return true;
  }

  void add_to_benefits(int vertex, int p) {
    if (p > 0) benefit(p - 1, p) = tournament(currentOrder[p - 1], vertex);

    for (int i = 0; i < p - 1; i++)
      benefit(i, p) = benefit(i, p - 1) + tournament(currentOrder[i], vertex);

    for (int i = 0; i < N; i++) outBegSecBenefit[i] -= tournament(i, vertex);
  }

  void delete_from_benefits(int vertex, int /*p*/) {
    for (int j = 0; j < N; j++) outBegSecBenefit[j] += tournament(j, vertex);
  }

  void init_benefits() {
    for (int i = 0; i < N; i++) {
      outBegSecBenefit[i] = 0;
      for (int j = 0; j < N; j++) outBegSecBenefit[i] += tournament(i, j);
    }
  }

  void branch_and_bound(int p, double begSecValue) {
    const auto tElapsed = timer.seconds();
    if (tElapsed > tMaxWalltime) {
      terminationReason = TerminationReason::MAX_WALL_TIME;
      return;  // max wall-time reached, clean-up
    }

    if (p == N - 1) {
      if (begSecValue < bound) {
        bound = begSecValue;
        std::copy_n(currentOrder.data(), N, bestOrder.data());
      }
      return;
    }

    const auto existResult = exist(p, begSecValue); /* BegSec Test */
    if (existResult == TreeTestResult::OPTIMAL_EARLY_TERMINATE) return;

    const auto q = N - p;
    const auto V = descent(q, currentOrder.data() + p);
    if (begSecValue + V < bound) {
      bound = begSecValue + V;
      std::copy_n(currentOrder.data(), N, bestOrder.data());
    }

    constexpr bool performRelaxation = true;
    if ((performRelaxation) && (q > 2)) { /*relax test */
      if ((existResult == TreeTestResult::RELAXATION_TEST) &&
          (bound - begSecValue <= lower_bound())) {
        return;
      }

      if (relax(currentOrder.data() + p, q, bound - begSecValue, V)) {
        return;
      }
    }

    auto &toExam_p = toExam[p];
    std::copy_n(currentOrder.data() + p, q, toExam_p.data());
    for (int i = 0; i < q; i++) {
      const auto vertex = toExam_p[i];
      if (!is_valid(vertex, p)) continue; /* ham, Smoves, OSmoves, lex tests */
      auto newBegSecValue = begSecValue;
      for (int j = p; j < N; j++) {
        if (tournament(vertex, currentOrder[j]) < 0)
          newBegSecValue += tournament(currentOrder[j], vertex);
      }
      exchange(vertex, p);
      add_to_benefits(vertex, p);
      isOutBegSec[vertex] = false;
      branch_and_bound(p + 1, newBegSecValue);
      isOutBegSec[vertex] = true;
      delete_from_benefits(vertex, p);
    }
  }

  void init_branch_and_bound() {
    toExam.resize(N);
    for (auto &vec : toExam) {
      vec.resize(N);
    }

    init_benefits();
    ro2 = ro1 = 2;
    head      = init_head_node();
  }

  void exchange(int vertex, int index) {
    int i = index;

    while (currentOrder[i] != vertex) i++;
    currentOrder[i]     = currentOrder[index];
    currentOrder[index] = vertex;
  }

  enum class TreeTestResult {
    NONE,
    OPTIMAL_EARLY_TERMINATE,
    RELAXATION_TEST,
  };

  TreeTestResult handle_leaf_node(std::shared_ptr<Node> &current, double begSecValue) {
    if (current->value < 0) {
      current->value = begSecValue;
      currentNode    = current;
      return TreeTestResult::NONE;
    }

    if (current->value <= begSecValue) {
      return TreeTestResult::OPTIMAL_EARLY_TERMINATE;
    } else {
      current->value = begSecValue;
      currentNode    = current;
      if (current->LowerBound >= 0) {
        return TreeTestResult::RELAXATION_TEST;
      } else {
        return TreeTestResult::NONE;
      }
    }

    return TreeTestResult::NONE;
  }

  TreeTestResult exist(int p, double BegSecValue) {
    int j         = -1;
    bool IsInTree = true;
    auto current  = head;

    for (int i = 0; i < p; i++) {
      while (isOutBegSec[++j])
        ;

      if (IsInTree) {
        process_current_node(current, j, IsInTree);

        if (i == p - 1) {
          return handle_leaf_node(current, BegSecValue);
        }
      } else {
        current = create(current, false, j);
        if (i == p - 1) current->value = BegSecValue;
      }
    }

    currentNode = current;
    return TreeTestResult::NONE;
  }

  void set_lower_bound(double value) {
    if (value > currentNode->LowerBound) currentNode->LowerBound = value;
  }

  [[nodiscard]] double lower_bound() const { return currentNode->LowerBound; }
};

long long unsigned int factorial(unsigned n) {
  long f = 1;
  for (unsigned i = 1; i <= n; ++i) f *= i;
  return f;
}

std::string gen_candidate_string(const std::vector<int> &solution) {
  std::ostringstream ss;
  bool first = true;
  ss << "(";
  std::for_each(solution.begin(), solution.end(), [&](const auto value) {
    if (!first) ss << ", ";
    ss << value;
    first = false;
  });
  ss << ")";
  return ss.str();
}

}  // namespace

std::pair<std::vector<int>, double> compute_min_ordering(
    const Teuchos::RCP<const Teuchos::Comm<int>> &communicator, double tMaxWalltime,
    const BlockNormsViewType &blockNorms, bool upperTriangular,
    std::function<double(const std::vector<int> &)> objectiveFunction,
    Teuchos::RCP<Teuchos::FancyOStream> out) {
  const auto tournament = convert_to_normal_form(blockNorms);

  Kokkos::Timer timer;
  LinearOrderingSolver solver{tournament, tMaxWalltime, communicator};
  solver.solve();
  const auto &bestOrder        = solver.order();
  const auto terminationReason = solver.termination_reason();
  const auto tElapsed          = timer.seconds();

  const auto N     = tournament.extent(0);
  const auto N_fac = factorial(N);

  std::vector<int> order(N);
  for (auto i = 0U; i < N; i++) {
    order[i] = bestOrder[i];
  }
  if (!upperTriangular) {
    std::reverse(order.begin(), order.end());
  }

  auto obj_f = std::sqrt(objectiveFunction(order));
  if (out) {
    *out << "Branch and bound terminate due to " << to_string(terminationReason) << " : {\n";
    *out << "  N               : " << N << "\n";
    *out << "  N!              : " << N_fac << "\n";
    *out << "  Obj. func       : " << obj_f << "\n";
    *out << "  Minimizer       : " << gen_candidate_string(order) << "\n";
    *out << "  tElapsed        : " << tElapsed << "\n";
    *out << "}\n";
  }

  return std::make_pair(order, obj_f);
}

}  // namespace Teko::LinearOrdering
