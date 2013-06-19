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
        int pID,
        int dim;
        int numCorners;
        scalar_t **corners;
        scalar_t *lmins, *gmins;
        scalar_t *lmaxs, *gmaxs;
        scalar_t maxScalar;
        coordinatePartBox(int pid, int dim_, scalar_t *lMins, scalar_t *gMins,
                                    scalar_t *lMaxs, scalar_t *gMaxs):
            pID(pid),
            dim(dim_),
            numCorners(int(pow(2, dim_))),
            corners(0),
            lmins(lMins), gmins(gMins), lmaxs(lMaxs), gmaxs(gMaxs),
            maxScalar (std::numeric_limits<scalar_t>::max()){
            this->corners = new scalar_t *[dim];
            for (int i = 0; i < dim; ++i){
                this->corners[i] = new scalar_t[this->numCorners];
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

        ~coordinatePartBox(){
            
            for (int i = 0; i < dim; ++i){
                delete [] this->corners[i];
            }
            delete [] this->corners[i];
            //delete []this->mins;
            //delete [] this->maxs;
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
        //expensive to do all parts one by one.
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

template <typename Adapter, typename partId_t>
class CoordinateCommGraph{
  private:
  
  const Environment *env;
  const Teuchos::Comm<int> *comm;
  const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords;
  const Zoltan2::PartitioningSolution<Adapter> *soln;
  vector<coordinatePartBox> cpb;
  int coordDim;
  partId_t numParts;
` 
  
  public:

	CoordinateCommGraph(
      const Environment *env_,
      const Teuchos::Comm<int> *comm_,
      const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords_,
      const Zoltan2::PartitioningSolution<Adapter> *soln_
	):env(env_), comm(comm_),coords(coords_), soln(soln_){
    coordDim = coords_->getCoordinateDim();
    numParts = this->soln->getActualGlobalNumberOfParts();

    this->create_part_boxes();
    this->hash_part_boxes();
    this->find_neighbors();
  }

  void create_part_boxes(){

    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;

    typedef typename Adapter::scalar_t scalar_t;

    scalar_t *lmins = new scalar_t [numParts * coordDim];
    scalar_t *gmins = new scalar_t [numParts * coordDim];
    scalar_t *lmaxs = new scalar_t [numParts * coordDim];
    scalar_t *gmaxs = new scalar_t [numParts * coordDim];

    for(partId_t i = 0; i < numParts; ++i){
      coordinatePartBox tmp(i, this->coordDim,
                            lmins + i * coordDim,
                            gmins + i * coordDim,
                            lmaxs + i * coordDim,
                            gmaxs + i * coordDim
                          );
      cpb.push_back(tmp);
    }
    
    typedef StridedData<lno_t, scalar_t> input_t;
    ArrayView<const gno_t> gnos;
    ArrayView<input_t>     xyz;
    ArrayView<input_t>     wgts;
    coords->getCoordinates(gnos, xyz, wgts);

    //local and global num coordinates.
    lno_t numLocalCoords = coords->getLocalNumCoordinates();

    scalar_t **pqJagged_coordinates = allocMemory<scalar_t *>(coordDim);

    for (int dim=0; dim < coordDim; dim++){
        ArrayRCP<const scalar_t> ar;
        xyz[dim].getInputArray(ar);
        //pqJagged coordinate values assignment
        pqJagged_coordinates[dim] =  (scalar_t *)ar.getRawPtr();
    }

    partId_t *sol_part = soln->getPartList();
    for(lno_t i = 0; i < numLocalCoords; ++i){
      partId_t p = sol_part[i];
      cpb[p].updateMinMax(pqJagged_coordinates, i);
    }
    delete []pqJagged_coordinates;
    

    reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MIN,
        dim * numParts, lmins, gmins
        );
    reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_MAX,
        dim * numParts, lmaxs, gmaxs
        );
  }
  
  void hash_part_boxes (){
    partId_t pSingleDim = pow(double(numParts), double(1.0 / coorDim));
    vector < vector <partId_t> > hash(pSingleDim * coordDim);
    scalar_t *allMins = new scalar_t [coordDim];
    scalar_t *allMaxs = new scalar_t [coordDim];
    
    for (int j = 0; j < coordDim; ++j){ 
      allMins[j] = cpb[0].gmins[j];
      allMaxs[j] = cpb[0].gmaxs[j];

    }

    for (partId_t i = 1; i < numParts; ++i){
      for (int j = 0; j < coordDim; ++j){
        scalar_t minC = cpb[i].gmins[i];
        scalar_t maxC = cpb[i].gmaxs[i];
        if (minC < allMins[j]) allMins[j] = minC;
        if (maxC > allMaxs[j]) allMaxs[j] = maxC;
      }
    }
    scalar_t *slices = new scalar_t [coordDim];
    for (int j = 0; j < coordDim; ++j){
      slices[j] = (allMaxs[j] - allMins[j]) / pSingleDim;
    }
    
    scalar_t *minIndices = new scalar_t [coordDim];
    scalar_t *maxIndices = new scalar_t [coordDim];
    
    for (partId_t i = 0; i < numParts; ++i){
      for (int j = 0; j < coordDim; ++j){
        scalar_t minC = cpb[i].gmins[i];
        scalar_t maxC = cpb[i].gmaxs[i];
        minIndices[i] = (minC - allMins[j]) / slice;
        maxIndices[i]  = (maxC - allMins[j]) / slice;
      }
      partId_t hashToAdd = 0;
      /*
      for (int j = 0; j < coordDim; ++j){
        for (partId_t i = ; i < )
      }
*/

    }

    
  }
  

};

} // namespace Zoltan2
