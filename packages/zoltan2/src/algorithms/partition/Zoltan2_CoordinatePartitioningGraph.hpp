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

#ifndef _ZOLTAN2_COORDCOMMGRAPH_HPP_
#define _ZOLTAN2_COORDCOMMGRAPH_HPP_


#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <set>
#include <fstream>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_RCPDecl.hpp"

namespace Zoltan2{


#define Z2_ABS(x) ((x) >= 0 ? (x) : -(x))

/*! \brief coordinateModelPartBox Class,
 * represents the boundaries of the box which is a result of a geometric partitioning algorithm.
 */
template <typename scalar_t,typename part_t>
class coordinateModelPartBox{

        part_t pID; //part Id
        int dim;    //dimension of the box
        scalar_t *lmins;    //minimum boundaries of the box along all dimensions.
        scalar_t *lmaxs;    //maximum boundaries of the box along all dimensions.
        scalar_t maxScalar;
        scalar_t _EPSILON;

        //to calculate the neighbors of the box and avoid the p^2 comparisons,
        //we use hashing. A box can be put into multiple hash buckets.
        //the following 2 variable holds the minimum and maximum of the
        //hash values along all dimensions.
        part_t *minHashIndices;
        part_t *maxHashIndices;

        //result hash bucket indices.
        std::vector <part_t> *gridIndices;
        //neighbors of the box.
        std::set <part_t> neighbors;
public:
        /*! \brief Constructor
         */
        coordinateModelPartBox(part_t pid, int dim_):
            pID(pid),
            dim(dim_),
            lmins(0), lmaxs(0),
            maxScalar (std::numeric_limits<scalar_t>::max()),
            _EPSILON(std::numeric_limits<scalar_t>::epsilon()),
            minHashIndices(0),
            maxHashIndices(0),
            gridIndices(0), neighbors()
        {
            lmins = new scalar_t [dim];
            lmaxs = new scalar_t [dim];

            minHashIndices = new part_t [dim];
            maxHashIndices = new part_t [dim];
            gridIndices = new std::vector <part_t> ();
            for (int i = 0; i < dim; ++i){
                lmins[i] = -this->maxScalar;
                lmaxs[i] = this->maxScalar;
            }
        }
        /*! \brief  Constructor
         * deep copy of the maximum and minimum boundaries.
         */
        coordinateModelPartBox(part_t pid, int dim_, scalar_t *lmi, scalar_t *lma):
            pID(pid),
            dim(dim_),
            lmins(0), lmaxs(0),
            maxScalar (std::numeric_limits<scalar_t>::max()),
            _EPSILON(std::numeric_limits<scalar_t>::epsilon()),
            minHashIndices(0),
            maxHashIndices(0),
            gridIndices(0), neighbors()
        {
            lmins = new scalar_t [dim];
            lmaxs = new scalar_t [dim];
            minHashIndices = new part_t [dim];
            maxHashIndices = new part_t [dim];
            gridIndices = new std::vector <part_t> ();
            for (int i = 0; i < dim; ++i){
                lmins[i] = lmi[i];
                lmaxs[i] = lma[i];
            }
        }


        /*! \brief  Copy Constructor
         * deep copy of the maximum and minimum boundaries.
         */
        coordinateModelPartBox(const coordinateModelPartBox <scalar_t, part_t> &other):
            pID(other.getpId()),
            dim(other.getDim()),
            lmins(0), lmaxs(0),
            maxScalar (std::numeric_limits<scalar_t>::max()),
            _EPSILON(std::numeric_limits<scalar_t>::epsilon()),
            minHashIndices(0),
            maxHashIndices(0),
            gridIndices(0), neighbors()
        {

            lmins = new scalar_t [dim];
            lmaxs = new scalar_t [dim];
            minHashIndices = new part_t [dim];
            maxHashIndices = new part_t [dim];
            gridIndices = new std::vector <part_t> ();
            scalar_t *othermins = other.getlmins();
            scalar_t *othermaxs = other.getlmaxs();
            for (int i = 0; i < dim; ++i){
                lmins[i] = othermins[i];
                lmaxs[i] = othermaxs[i];
            }
        }
        /*! \brief  Destructor
         */
        ~coordinateModelPartBox(){
            delete []this->lmins;
            delete [] this->lmaxs;
            delete []this->minHashIndices;
            delete [] this->maxHashIndices;
            delete gridIndices;
        }

        /*! \brief  function to set the part id
         */
        void setpId(part_t pid){
            this->pID = pid;
        }
        /*! \brief  function to get the part id
         */
        part_t getpId() const{
            return this->pID;
        }


        /*! \brief  function to set the dimension
         */
        int getDim()const{
            return this->dim;
        }
        /*! \brief  function to get minimum values along all dimensions
         */
        scalar_t * getlmins()const{
            return this->lmins;
        }
        /*! \brief  function to get maximum values along all dimensions
         */
        scalar_t * getlmaxs()const{
            return this->lmaxs;
        }
        /*! \brief  compute the centroid of the box 
         */
        void computeCentroid(scalar_t *centroid)const {
            for (int i = 0; i < this->dim; i++)
                centroid[i] = 0.5 * (this->lmaxs[i] + this->lmins[i]);
        }

        /*! \brief  function to get the indices of the buckets
         * that the part is inserted to
         */
        std::vector <part_t> * getGridIndices (){
            return this->gridIndices;
        }

        /*! \brief  function to get the indices of the neighboring parts.
         */
        std::set<part_t> *getNeighbors(){
            return &(this->neighbors);
        }

        /*! \brief function to test whether a point is in the box
         */
        bool pointInBox(int pointdim, scalar_t *point) {
          if (pointdim != this->dim) 
            throw std::logic_error("dim of point must match dim of box");
          for (int i = 0; i < pointdim; i++) {
            if (point[i] < this->lmins[i]) return false;
            if (point[i] > this->lmaxs[i]) return false;
          }
          return true;
        }

        /*! \brief  function to check if two boxes are neighbors.
         */
        bool isNeighborWith(const coordinateModelPartBox <scalar_t, part_t> &other){


            scalar_t *omins = other.getlmins();
            scalar_t *omaxs = other.getlmaxs();

            int equality = 0;
            for (int i = 0; i < dim; ++i){

                if (omins[i] - this->lmaxs[i] > _EPSILON  || this->lmins[i] - omaxs[i] > _EPSILON ){
                    return false;
                }
                else if (Z2_ABS(omins[i] - this->lmaxs[i]) < _EPSILON  || Z2_ABS(this->lmins[i] - omaxs[i]) < _EPSILON ){
                    if (++equality > 1){
                        return false;
                    }
                }
            }
            if (equality == 1) {
                return true;
            }
            else {
                std::cout << "something is wrong: equality:" << equality << std::endl;
                return false;
            }
        }


        /*! \brief  function to add a new neighbor to the neighbor list.
         */
        void addNeighbor(part_t nIndex){
            neighbors.insert(nIndex);
        }
        /*! \brief  function to check if a given part is already in the neighbor list.
         */
        bool isAlreadyNeighbor(part_t nIndex){

            if (neighbors.end() != neighbors.find(nIndex)){
                return true;
            }
            return false;

        }


        /*! \brief  function to obtain the min and max hash values along all dimensions.
         */
        void setMinMaxHashIndices (
                scalar_t *minMaxBoundaries,
                scalar_t *sliceSizes,
                part_t numSlicePerDim
                ){
            for (int j = 0; j < dim; ++j){
                scalar_t distance = (lmins[j] - minMaxBoundaries[j]);
                part_t minInd = 0;
                if (distance > _EPSILON && sliceSizes[j] > _EPSILON){
                    minInd = static_cast<part_t>(floor((lmins[j] - minMaxBoundaries[j])/ sliceSizes[j]));
                }

                if(minInd >= numSlicePerDim){
                    minInd = numSlicePerDim - 1;
                }
                part_t maxInd = 0;
                distance = (lmaxs[j] - minMaxBoundaries[j]);
                if (distance > _EPSILON && sliceSizes[j] > _EPSILON){
                    maxInd = static_cast<part_t>(ceil((lmaxs[j] - minMaxBoundaries[j])/ sliceSizes[j]));
                }
                if(maxInd >= numSlicePerDim){
                    maxInd = numSlicePerDim - 1;
                }

                //cout << "j:" << j << " lmins:" << lmins[j] << " lmaxs:" << lmaxs[j] << endl;
                //cout << "j:" << j << " min:" << minInd << " max:" << maxInd << endl;
                minHashIndices[j] = minInd;
                maxHashIndices[j] = maxInd;
            }
            std::vector <part_t> *in = new std::vector <part_t> ();
            in->push_back(0);
            std::vector <part_t> *out = new std::vector <part_t> ();

            for (int j = 0; j < dim; ++j){

                part_t minInd = minHashIndices[j];
                part_t maxInd = maxHashIndices[j];


                part_t pScale = part_t(pow (float(numSlicePerDim), int(dim - j -1)));

                part_t inSize = in->size();

                for (part_t k = minInd; k <= maxInd; ++k){
                    for (part_t i = 0; i < inSize; ++i){
                        out->push_back((*in)[i] + k * pScale);
                    }
                }
                in->clear();
                std::vector <part_t> *tmp = in;
                in= out;
                out= tmp;
            }

            std::vector <part_t> *tmp = in;
            in = gridIndices;
            gridIndices = tmp;


            delete in;
            delete out;
        }

        /*! \brief  function to print the boundaries.
        */
        void print(){
            for(int i = 0; i < this->dim; ++i){
                std::cout << "\tbox:" << this->pID << " dim:" << i << " min:" << lmins[i] << " max:" << lmaxs[i] << std::endl;
            }
        }

        /*! \brief  function to update the boundary of the box.
        */
        void updateMinMax (scalar_t newBoundary, int isMax, int dimInd){
            if (isMax){
                lmaxs[dimInd] = newBoundary;
            }
            else {
                lmins[dimInd] = newBoundary;
            }
        }

        template <typename tt>
        std::string toString(tt obj){
            std::stringstream ss (std::stringstream::in |std::stringstream::out);
            ss << obj;
            std::string tmp = "";
            ss >> tmp;
            return tmp;
        }


        /*! \brief  function for visualization.
        */
        void writeGnuPlot(std::ofstream &file,std::ofstream &mm){
            int numCorners = (int(1)<<dim);
            scalar_t *corner1 = new scalar_t [dim];
            scalar_t *corner2 = new scalar_t [dim];

            for (int i = 0; i < dim; ++i){
                /*
                if (-maxScalar == lmins[i]){
                    if (lmaxs[i] > 0){
                        lmins[i] = lmaxs[i] / 2;
                    }
                    else{
                        lmins[i] = lmaxs[i] * 2;
                    }
                }
                */
                //std::cout << lmins[i] << " ";
                mm << lmins[i] << " ";
            }
            //std::cout <<  std::endl;
            mm << std::endl;
            for (int i = 0; i < dim; ++i){

                /*
                if (maxScalar == lmaxs[i]){
                    if (lmins[i] < 0){
                        lmaxs[i] = lmins[i] / 2;
                    }
                    else{
                        lmaxs[i] = lmins[i] * 2;
                    }
                }
                */

                //std::cout << lmaxs[i] << " ";
                mm << lmaxs[i] << " ";
            }
            //std::cout <<  std::endl;
            mm << std::endl;

            for (int j = 0; j < numCorners; ++j){
                std::vector <int> neighborCorners;
                for (int i = 0; i < dim; ++i){
                    if(int(j & (int(1)<<i)) == 0){
                        corner1[i] = lmins[i];
                    }
                    else {
                        corner1[i] = lmaxs[i];
                    }
                    if (j % (int(1)<<(i + 1)) >= (int(1)<<i)){
                        int c1 = j - (int(1)<<i);

                        if (c1 > 0) {
                            neighborCorners.push_back(c1);
                        }
                    }
                    else {

                        int c1 = j + (int(1)<<i);
                        if (c1 < (int(1) << dim)) {
                            neighborCorners.push_back(c1);
                        }
                    }
                }
                //std::cout << "me:" << j << " nc:" << int (neighborCorners.size()) << std::endl;
                for (int m = 0; m < int (neighborCorners.size()); ++m){

                    int n = neighborCorners[m];
                    //std::cout << "me:" << j << " n:" << n << std::endl;
                    for (int i = 0; i < dim; ++i){
                        if(int(n & (int(1)<<i)) == 0){
                            corner2[i] = lmins[i];
                        }
                        else {
                            corner2[i] = lmaxs[i];
                        }
                    }

                    std::string arrowline = "set arrow from ";
                    for (int i = 0; i < dim - 1; ++i){
                        arrowline += toString<scalar_t>(corner1[i]) + ",";
                    }
                    arrowline += toString<scalar_t>(corner1[dim -1]) + " to ";

                    for (int i = 0; i < dim - 1; ++i){
                        arrowline += toString<scalar_t>(corner2[i]) + ",";
                    }
                    arrowline += toString<scalar_t>(corner2[dim -1]) + " nohead\n";

                    file << arrowline;
                }
            }
            delete []corner1;
            delete []corner2;
        }


};


/*! \brief GridHash Class,
 * Hashing Class for part boxes
 */
template <typename scalar_t, typename part_t>
class GridHash{
private:

    RCP < std::vector <Zoltan2::coordinateModelPartBox <scalar_t, part_t> > > pBoxes;

    //minimum of the maximum box boundaries
    scalar_t *minMaxBoundaries;
    //maximum of the minimum box boundaries
    scalar_t *maxMinBoundaries;
    //the size of each slice along dimensions
    scalar_t *sliceSizes;
    part_t nTasks;
    int dim;
    //the number of slices per dimension
    part_t numSlicePerDim;
    //the number of grids - buckets
    part_t numGrids;
    //hash vector
    std::vector <std::vector <part_t>  > grids;
    //result communication graph.
    ArrayRCP <part_t> comXAdj;
    ArrayRCP <part_t> comAdj;
public:

    /*! \brief GridHash Class,
     * Constructor
     */
    GridHash(RCP < std::vector <Zoltan2::coordinateModelPartBox <scalar_t, part_t> > > pBoxes_,
            part_t ntasks_, int dim_):
        pBoxes(pBoxes_),
        minMaxBoundaries(0),
        maxMinBoundaries(0), sliceSizes(0),
        nTasks(ntasks_),
        dim(dim_),
        numSlicePerDim(part_t(pow(double(ntasks_), 1.0 / dim))),
        numGrids(0),
        grids(),
        comXAdj(), comAdj()
    {

        minMaxBoundaries = new scalar_t[dim];
        maxMinBoundaries = new scalar_t[dim];
        sliceSizes = new scalar_t[dim];
        //calculate the number of slices in each dimension.
        numSlicePerDim /= 2;
        if (numSlicePerDim == 0) numSlicePerDim = 1;

        numGrids = part_t(pow(float(numSlicePerDim), int(dim)));

        //allocate memory for buckets.
        std::vector <std::vector <part_t>  > grids_ (numGrids);
        this->grids = grids_;
        //get the boundaries of buckets.
        this->getMinMaxBoundaries();
        //insert boxes to buckets
        this->insertToHash();
        //calculate the neighbors for each bucket.
        part_t nCount = this->calculateNeighbors();

        //allocate memory for communication graph
        ArrayRCP <part_t> tmpComXadj(ntasks_);
        ArrayRCP <part_t> tmpComAdj(nCount);
        comXAdj = tmpComXadj;
        comAdj = tmpComAdj;
        //fill communication graph
        this->fillAdjArrays();
    }


    /*! \brief GridHash Class,
     * Destructor
     */
    ~GridHash(){
        delete []minMaxBoundaries;
        delete []maxMinBoundaries;
        delete []sliceSizes;
    }

    /*! \brief GridHash Class,
     * Function to fill adj arrays.
     */
    void fillAdjArrays(){

        part_t adjIndex = 0;

        for(part_t i = 0; i < this->nTasks; ++i){
            std::set<part_t> *neigbors = (*pBoxes)[i].getNeighbors();

            part_t s = neigbors->size();

            comXAdj[i] = s;
            if (i > 0){
                comXAdj[i] += comXAdj[i - 1];
            }
            typedef typename std::set<part_t> mySet;
            typedef typename mySet::iterator myIT;
            myIT it;
            for (it=neigbors->begin(); it!=neigbors->end(); ++it)

                comAdj[adjIndex++] = *it;
            //TODO not needed anymore.
            neigbors->clear();
        }
    }



    /*! \brief GridHash Class,
     * returns the adj arrays.
     */
    void getAdjArrays(
            ArrayRCP <part_t> &comXAdj_,
            ArrayRCP <part_t> &comAdj_){
        comXAdj_ = this->comXAdj;
        comAdj_ = this->comAdj;
    }

    /*! \brief GridHash Class,
     * For each box compares the adjacency against the boxes that are in the same buckets.
     */
    part_t calculateNeighbors(){
        part_t nCount = 0;
        for(part_t i = 0; i < this->nTasks; ++i){
            std::vector <part_t> *gridIndices =(*pBoxes)[i].getGridIndices();
            part_t gridCount = gridIndices->size();

            for (part_t j = 0; j < gridCount; ++j){
                part_t grid = (*gridIndices)[j];
                part_t boxCount = grids[grid].size();
                for (part_t k = 0; k < boxCount; ++k){
                    part_t boxIndex = grids[grid][k];
                    if (boxIndex > i){
                        if((!(*pBoxes)[i].isAlreadyNeighbor(boxIndex))&& (*pBoxes)[i].isNeighborWith((*pBoxes)[boxIndex])){
                            //cout << "i:" << i << " n:" << boxIndex << " are neighbors."<< endl;
                            (*pBoxes)[i].addNeighbor(boxIndex);
                            (*pBoxes)[boxIndex].addNeighbor(i);
                            nCount += 2;
                        }
                    }
                }
            }
        }

        return nCount;
    }

    /*! \brief GridHash Class,
     * For each box calculates the buckets which it should be inserted to.
     */
    void insertToHash(){

        //cout << "ntasks:" << this->nTasks << endl;
        for(part_t i = 0; i < this->nTasks; ++i){
            (*pBoxes)[i].setMinMaxHashIndices(minMaxBoundaries, sliceSizes, numSlicePerDim);

            std::vector <part_t> *gridIndices =(*pBoxes)[i].getGridIndices();

            part_t gridCount = gridIndices->size();
            //cout << "i:" << i << " gridsize:" << gridCount << endl;
            for (part_t j = 0; j < gridCount; ++j){
                part_t grid = (*gridIndices)[j];

                //cout << "i:" << i << " is being inserted to:" << grid << endl;
                (grids)[grid].push_back(i);
            }
        }


/*
        for(part_t i = 0; i < grids.size(); ++i){
            cout << "grid:" << i << " gridsuze:" << (grids)[i].size() << " elements:";
            for(part_t j = 0; j < (grids)[i].size(); ++j){
                cout <<(grids)[i][j] << " ";
            }
            cout << endl;

        }
*/
    }

    /*! \brief GridHash Class,
     * calculates the minimum of maximum box boundaries, and maxium of minimum box boundaries.
     */
    void getMinMaxBoundaries(){
        scalar_t *mins = (*pBoxes)[0].getlmins();
        scalar_t *maxs = (*pBoxes)[0].getlmaxs();

        for (int j = 0; j < dim; ++j){
            minMaxBoundaries[j] = maxs[j];
            maxMinBoundaries[j] = mins[j];
        }

        for (part_t i = 1; i < nTasks; ++i){

            mins = (*pBoxes)[i].getlmins();
            maxs = (*pBoxes)[i].getlmaxs();

            for (int j = 0; j < dim; ++j){

                if (minMaxBoundaries[j] > maxs[j]){
                    minMaxBoundaries[j] = maxs[j];
                }
                if (maxMinBoundaries[j] < mins[j]){
                    maxMinBoundaries[j] = mins[j];
                }
            }
        }


        for (int j = 0; j < dim; ++j){
            sliceSizes[j] = (maxMinBoundaries[j] - minMaxBoundaries[j]) / numSlicePerDim;
            if (sliceSizes[j] < 0) sliceSizes[j] = 0;
                /*
            cout << "dim:" << j <<
                    " minMax:" <<  minMaxBoundaries[j] <<
                    " maxMin:" << maxMinBoundaries[j] <<
                    " sliceSizes:" << sliceSizes[j] << endl;
                    */
        }
    }
};
/*
template <typename scalar_t,typename part_t>
class coordinatePartBox{
public:
        part_t pID;
        int dim;
        int numCorners;
        scalar_t **corners;
        scalar_t *lmins, *gmins;
        scalar_t *lmaxs, *gmaxs;
        scalar_t maxScalar;
        std::vector <part_t> hash_indices;
        coordinatePartBox(part_t pid, int dim_, scalar_t *lMins, scalar_t *gMins,
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
                lmins[i] = this->maxScalar;
                lmaxs[i] = -this->maxScalar;
            }


            for (int j = 0; j < this->numCorners; ++j){
                for (int i = 0; i < dim; ++i){
                    std::cout << "j:" << j << " i:" << i << " 2^i:" << pow(2,i) << " and:" << int(j & int(pow(2,i))) << std::endl;
                    if(int(j & int(pow(2,i))) == 0){
                        corners[i][j] = gmins[i];
                    }
                    else {
                        corners[i][j] = gmaxs[i];
                    }

                }
            }
        }

};

template <typename Adapter, typename part_t>
class CoordinateCommGraph{
private:

    typedef typename Adapter::lno_t lno_t;
    typedef typename Adapter::gno_t gno_t;
    typedef typename Adapter::scalar_t scalar_t;

    const Environment *env;
    const Teuchos::Comm<int> *comm;
    const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords;
    const Zoltan2::PartitioningSolution<Adapter> *soln;
    std::vector<coordinatePartBox, part_t> cpb;
    int coordDim;
    part_t numParts;


public:

    CoordinateCommGraph(
            const Environment *env_,
            const Teuchos::Comm<int> *comm_,
            const Zoltan2::CoordinateModel<typename Adapter::base_adapter_t> *coords_,
            const Zoltan2::PartitioningSolution<Adapter> *soln_
    ):
        env(env_),
        comm(comm_),
        coords(coords_),
        soln(soln_),
        coordDim (coords_->getCoordinateDim()),
        numParts (this->soln->getActualGlobalNumberOfParts())
        {
        this->create_part_boxes();
        this->hash_part_boxes();
        this->find_neighbors();
    }

    void create_part_boxes(){


        size_t allocSize = numParts * coordDim;
        scalar_t *lmins = new scalar_t [allocSize];
        scalar_t *gmins = new scalar_t [allocSize];
        scalar_t *lmaxs = new scalar_t [allocSize];
        scalar_t *gmaxs = new scalar_t [allocSize];

        for(part_t i = 0; i < numParts; ++i){
            coordinatePartBox tmp(
                    i,
                    this->coordDim,
                    lmins + i * coordDim,
                    gmins + i * coordDim,
                    lmaxs + i * coordDim,
                    gmaxs + i * coordDim
            );
            cpb.push_back(tmp);
        }

        typedef StridedData<lno_t, scalar_t> input_t;
        Teuchos::ArrayView<const gno_t> gnos;
        Teuchos::ArrayView<input_t>     xyz;
        Teuchos::ArrayView<input_t>     wgts;
        coords->getCoordinates(gnos, xyz, wgts);

        //local and global num coordinates.
        lno_t numLocalCoords = coords->getLocalNumCoordinates();

        scalar_t **pqJagged_coordinates = new scalar_t *[coordDim];

        for (int dim=0; dim < coordDim; dim++){
            Teuchos::ArrayRCP<const scalar_t> ar;
            xyz[dim].getInputArray(ar);
            //pqJagged coordinate values assignment
            pqJagged_coordinates[dim] =  (scalar_t *)ar.getRawPtr();
        }

        part_t *sol_part = soln->getPartList();
        for(lno_t i = 0; i < numLocalCoords; ++i){
            part_t p = sol_part[i];
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
        part_t pSingleDim = pow(double(numParts), double(1.0 / coordDim));
        if (pSingleDim == 0) pSingleDim = 1;
        std::vector < std::vector <part_t> > hash
                (
                        part_t ( pow ( part_t (pSingleDim),
                                         part_t(coordDim)
                                       )
                                 )
                );

        //calculate the corners of the dataset.
        scalar_t *allMins = new scalar_t [coordDim];
        scalar_t *allMaxs = new scalar_t [coordDim];
        part_t *hash_scales= new scalar_t [coordDim];

        for (int j = 0; j < coordDim; ++j){
            allMins[j] = cpb[0].gmins[j];
            allMaxs[j] = cpb[0].gmaxs[j];
            hash_scales[j] = part_t ( pow ( part_t (pSingleDim), part_t(coordDim - j - 1)));
        }

        for (part_t i = 1; i < numParts; ++i){
            for (int j = 0; j < coordDim; ++j){
                scalar_t minC = cpb[i].gmins[i];
                scalar_t maxC = cpb[i].gmaxs[i];
                if (minC < allMins[j]) allMins[j] = minC;
                if (maxC > allMaxs[j]) allMaxs[j] = maxC;
            }
        }

        //get size of each hash for each dimension
        scalar_t *hash_slices_size = new scalar_t [coordDim];
        for (int j = 0; j < coordDim; ++j){
            hash_slices_size[j] = (allMaxs[j] - allMins[j]) / pSingleDim;

        }

        delete []allMaxs;
        delete []allMins;



        std::vector <part_t> *hashIndices = new std::vector <part_t>();
        std::vector <part_t> *resultHashIndices = new std::vector <part_t>();
        std::vector <part_t> *tmp_swap;
        for (part_t i = 0; i < numParts; ++i){
            hashIndices->clear();
            resultHashIndices->clear();

            hashIndices->push_back(0);

            for (int j = 0; j < coordDim; ++j){

                scalar_t minC = cpb[i].gmins[i];
                scalar_t maxC = cpb[i].gmaxs[i];
                part_t minHashIndex = part_t ((minC - allMins[j]) / hash_slices_size[j]);
                part_t maxHashIndex  = part_t ((maxC - allMins[j]) / hash_slices_size[j]);

                part_t hashIndexSize = hashIndices->size();

                for (part_t k = minHashIndex; k <= maxHashIndex; ++k ){

                    for (part_t i = 0; i < hashIndexSize; ++i){
                        resultHashIndices->push_back(hashIndices[i] + k *  hash_scales[j]);
                    }
                }
                tmp_swap = hashIndices;
                hashIndices = resultHashIndices;
                resultHashIndices = tmp_swap;
            }

            part_t hashIndexSize = hashIndices->size();
            for (part_t j = 0; j < hashIndexSize; ++j){
                hash[(*hashIndices)[j]].push_back(i);
            }
            cpb[i].hash_indices = (*hashIndices);
        }
        delete hashIndices;
        delete resultHashIndices;
    }

    void find_neighbors(){

    }


};

*/
} // namespace Zoltan2

#endif
