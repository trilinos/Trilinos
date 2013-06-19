// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER


#include <cmath>
#include <limits>
#include <iostream>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Comm.hpp"

namespace Zoltan2{
template <typename scalar_t>
class coordinatePartBox{
public:
        int dim;
        int numCorners;
        scalar_t **corners;
        scalar_t *mins;
        scalar_t *maxs;
        scalar_t maxScalar;
        coordinatePartBox(int dim_):
            dim(dim_),
            numCorners(int(pow(2, dim_))),
            corners(0),
            maxScalar (std::numeric_limits<scalar_t>::max()){
            this->corners = new scalar_t *[dim];
            for (int i = 0; i < dim; ++i){
                this->corners[i] = new scalar_t[this->numCorners];
            }
            this->mins = new scalar_t [this->dim];
            this->maxs = new scalar_t [this->dim];
            for (int i = 0; i < dim; ++i){
                this->mins[i] = this->maxScalar;
                this->maxs[i] = -this->maxScalar;
            }

            /*
            for (int j = 0; j < this->numCorners; ++j){
                for (int i = 0; i < dim; ++i){
                    if(j & pow(2,i) == 0){
                        corners[i][j] = this->maxScalar;
                    }
                    else {
                        corners[i][j] = -this->maxScalar;
                    }

                }
            }
            */

        }
        template <typename lno_t>
        void updateMinMax (scalar_t **coords, lno_t cIndex){
            for (int i = 0; i < dim; ++i){
                if (coords[i][cIndex] < mins[i]){
                    mins[i] = coords[i][cIndex];
                }
                if (coords[i][cIndex] > maxs[i]){
                    maxs[i] = coords[i][cIndex];
                }

            }
        }

        void getGlobalMinMax(Comm<int> *comm){

            scalar_t *global_mins = new scalar_t [this->dim];
            reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN,
                    dim, mins, global_mins
            );
            scalar_t *global_max = mins;
            this->mins = global_mins;
            reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MAX,
                    dim, maxs, global_max
            );

            delete [] maxs;
            maxs = global_max;
        }



        void updateCorners(){
            /*
            for (int i = 0; i < dim; ++i){
                std::cout << "i:" << i << " min:" << mins[i] << " max:" << maxs[i] << std::endl;
            }
            */
            for (int j = 0; j < this->numCorners; ++j){
                for (int i = 0; i < dim; ++i){
                    std::cout << "j:" << j << " i:" << i << " 2^i:" << pow(2,i) << " and:" << int(j & int(pow(2,i))) << std::endl;
                    if(int(j & int(pow(2,i))) == 0){
                        corners[i][j] = mins[i];
                    }
                    else {
                        corners[i][j] = maxs[i];
                    }

                }
            }
        }

};

template <typename Adapter>
class CoordinateCommGraph{
	CoordinateCommGraph(){
	}

};

} // namespace Zoltan2
