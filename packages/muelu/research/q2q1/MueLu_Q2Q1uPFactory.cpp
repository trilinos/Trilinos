// @HEADER
//
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_Q2Q1UPFACTORY_DECL_HPP
#define MUELU_Q2Q1UPFACTORY_DECL_HPP

#include <Teuchos_RCP.hpp>

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_Q2Q1uPFactory.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"

#include <algorithm>
#include <vector>
#include <sys/stat.h>

namespace MueLu {

  std::string i2s(int i) {
    std::ostringstream os;
    if (i < 10)    os << "0";
    if (i < 100)   os << "0";
    if (i < 1000)  os << "0";
    if (i < 10000) os << "0";
    os << i;
    return os.str();
  }

  // Sort a double array and move along list2 to match sorted array
  void Muelu_az_dsort2(std::vector<double>& dlist, std::vector<int>& list2) {
    int    l, r, j, i, flag;
    int    RR2;
    double dRR, dK;

    int N = dlist.size();
    if (N <= 1) return;

    l    = N / 2 + 1;
    r    = N - 1;
    l    = l - 1;
    dRR  = dlist[l - 1];
    dK   = dlist[l - 1];

    if (list2.size()) {
      RR2 = list2[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (dlist[j] > dlist[j - 1]) j = j + 1;

            if (dlist[j - 1] > dK) {
              dlist[i - 1] = dlist[j - 1];
              list2[i - 1] = list2[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }
        dlist[i - 1] = dRR;
        list2[i - 1] = RR2;

        if (l == 1) {
          dRR = dlist[r];
          RR2 = list2[r];
          dK = dlist[r];
          dlist[r] = dlist[0];
          list2[r] = list2[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          dRR = dlist[l - 1];
          RR2 = list2[l - 1];
          dK  = dlist[l - 1];
        }
      }
      dlist[0] = dRR;
      list2[0] = RR2;
    }
    else {
      while (r != 0) {
        j = l;
        flag = 1;
        while (flag == 1) {
          i = j;
          j = j + j;
          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (dlist[j] > dlist[j - 1]) j = j + 1;
            if (dlist[j - 1] > dK) {
              dlist[i - 1] = dlist[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }
        dlist[i - 1] = dRR;
        if (l == 1) {
          dRR = dlist[r];
          dK  = dlist[r];
          dlist[r] = dlist[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          dRR  = dlist[l - 1];
          dK   = dlist[l - 1];
        }
      }
      dlist[0] = dRR;
    }
  }

  /* ******************************************************************* */
  /* sort an array and move along list2 and/or list to match sorted array*/
  /* ------------------------------------------------------------------- */
  void Muelu_az_sort(int list[], int N, int list2[], double list3[]) {
    int    l, r, RR, K, j, i, flag;
    int    RR2;
    double RR3;

    if (N <= 1) return;

    l   = N / 2 + 1;
    r   = N - 1;
    l   = l - 1;
    RR  = list[l - 1];
    K   = list[l - 1];

    if ((list2 != NULL) && (list3 != NULL)) {
      RR2 = list2[l - 1];
      RR3 = list3[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list [i - 1] = list [j - 1];
              list2[i - 1] = list2[j - 1];
              list3[i - 1] = list3[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list [i - 1] = RR;
        list2[i - 1] = RR2;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list [r];
          RR2 = list2[r];
          RR3 = list3[r];

          K = list[r];
          list[r ] = list[0];
          list2[r] = list2[0];
          list3[r] = list3[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list [l - 1];
          RR2 = list2[l - 1];
          RR3 = list3[l - 1];
          K   = list [l - 1];
        }
      }

      list [0] = RR;
      list2[0] = RR2;
      list3[0] = RR3;
    }
    else if (list2 != NULL) {
      RR2 = list2[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list [i - 1] = list [j - 1];
              list2[i - 1] = list2[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list [i - 1] = RR;
        list2[i - 1] = RR2;

        if (l == 1) {
          RR  = list [r];
          RR2 = list2[r];

          K = list[r];
          list[r ] = list[0];
          list2[r] = list2[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list [l - 1];
          RR2 = list2[l - 1];
          K   = list [l - 1];
        }
      }

      list [0] = RR;
      list2[0] = RR2;
    }
    else if (list3 != NULL) {
      RR3 = list3[l - 1];
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list [i - 1] = list [j - 1];
              list3[i - 1] = list3[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list [i - 1] = RR;
        list3[i - 1] = RR3;

        if (l == 1) {
          RR  = list [r];
          RR3 = list3[r];

          K = list[r];
          list[r ] = list[0];
          list3[r] = list3[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list [l - 1];
          RR3 = list3[l - 1];
          K   = list [l - 1];
        }
      }

      list [0] = RR;
      list3[0] = RR3;

    }
    else {
      while (r != 0) {
        j = l;
        flag = 1;

        while (flag == 1) {
          i = j;
          j = j + j;

          if (j > r + 1)
            flag = 0;
          else {
            if (j < r + 1)
              if (list[j] > list[j - 1]) j = j + 1;

            if (list[j - 1] > K) {
              list[i - 1] = list[j - 1];
            }
            else {
              flag = 0;
            }
          }
        }

        list[i - 1] = RR;

        if (l == 1) {
          RR  = list [r];

          K = list[r];
          list[r ] = list[0];
          r = r - 1;
        }
        else {
          l   = l - 1;
          RR  = list[l - 1];
          K   = list[l - 1];
        }
      }

      list[0] = RR;
    }
  }

  // Merge two already sorted lists into one combined sorted list.
  // NOTE: lists are given as integer arrays. These integer arrays give
  // locations in CoordDist[] defining the list values. That the ith value
  // associated with the Candidates list is actually CoordDist[Candidates[i]].
  void MergeSort(std::vector<int>& oldCandidates, size_t numOldCandidates, const std::vector<int>& newCandidates, const std::vector<double>& coordDist, ArrayRCP<const size_t> ia) {
    size_t numNewCandidates = newCandidates.size();
    size_t numCandidates    = numOldCandidates + numNewCandidates;

    oldCandidates.resize(numCandidates);

    int i = numOldCandidates - 1;
    int j = numNewCandidates - 1;
    int k = numCandidates    - 1;
    while ((i >= 0) || (j >= 0)) {
      if      (i < 0) oldCandidates[k--] = newCandidates[j--];
      else if (j < 0) oldCandidates[k--] = oldCandidates[i--];
      else {
        int ii = oldCandidates[i];
        int jj = newCandidates[j];

        // Must match code above. There is something arbitrary
        // and crappy about the current weighting.

#ifdef optimal
        if (-coordDist[ii] - .01*(ia[ii+1]-ia[ii]) + 1.e-10*(ii+1) <
            -coordDist[jj] - .01*(ia[jj+1]-ia[jj]) + 1.e-10*(jj+1))
#else
        if (coordDist[ii] - .0*(ia[ii+1]-ia[ii]) + 1.e-3*(ii+1) <
            coordDist[jj] - .0*(ia[jj+1]-ia[jj]) + 1.e-3*(jj+1))
        // if (ii < jj)
#endif
          oldCandidates[k--] = newCandidates[j--];
        else
          oldCandidates[k--] = oldCandidates[i--];
      }
    }
  }

} // namespace MueLu

#endif // MUELU_Q2Q1UPFACTORY_DECL_HPP
