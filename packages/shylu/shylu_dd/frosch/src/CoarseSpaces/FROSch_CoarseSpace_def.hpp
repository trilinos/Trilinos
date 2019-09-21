//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_COARSESPACE_DEF_HPP
#define _FROSCH_COARSESPACE_DEF_HPP

#include <FROSch_CoarseSpace_decl.hpp>


namespace FROSch {
    
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    CoarseSpace<SC,LO,GO,NO>::CoarseSpace() :
    SerialRowMap_ (),
    UnassembledBasesMaps_ (0),
    UnassembledSubspaceBases_ (0),
    AssembledBasisMap_ (),
    AssembledBasis_ (),
    GlobalBasisMatrix_ ()
    {

    }

    // Will man Informationen über die Subspaces als strings reingeben?
    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::addSubspace(XMapPtr subspaceBasisMap,
                                              XMultiVectorPtr localSubspaceBasis)
    {
        FROSCH_ASSERT(!subspaceBasisMap.is_null(),"subspaceBasisMap.is_null()");
        if (!localSubspaceBasis.is_null()) {
            FROSCH_ASSERT(localSubspaceBasis->getNumVectors()==subspaceBasisMap->getNodeNumElements(),"localSubspaceBasis->getNumVectors()!=subspaceBasisMap->getNodeNumElements()");
            if (!SerialRowMap_.is_null()) {
                FROSCH_ASSERT(SerialRowMap_->isSameAs(*localSubspaceBasis->getMap()),"!UnassembledSubspaceBases_[0]->isSameAs(localSubspaceBasis->getMap())");
            } else {
                SerialRowMap_ = localSubspaceBasis->getMap();
            }
        } else {
            FROSCH_ASSERT(subspaceBasisMap->getNodeNumElements()==0,"subspaceBasisMap->getNodeNumElements()!=0");
        }
        UnassembledBasesMaps_.push_back(subspaceBasisMap);
        UnassembledSubspaceBases_.push_back(localSubspaceBasis);
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::assembleCoarseSpace()
    {
        FROSCH_ASSERT(UnassembledBasesMaps_.size()>0,"UnassembledBasesMaps_.size()==0");
        FROSCH_ASSERT(UnassembledSubspaceBases_.size()>0,"UnassembledSubspaceBases_.size()==0");

        UN itmp = 0;
        LOVecPtr2D partMappings;
        AssembledBasisMap_ = AssembleMaps(UnassembledBasesMaps_(),partMappings);
        if (!AssembledBasisMap_.is_null()&&!SerialRowMap_.is_null()) {
            if (AssembledBasisMap_->getGlobalNumElements()>0) { // AH 02/12/2019: Is this the right condition? Seems to work for now...
                AssembledBasis_ = MultiVectorFactory<SC,LO,GO,NO >::Build(SerialRowMap_,AssembledBasisMap_->getNodeNumElements());
                for (UN i=0; i<UnassembledBasesMaps_.size(); i++) {
                    for (UN j=0; j<UnassembledBasesMaps_[i]->getNodeNumElements(); j++) {
                        AssembledBasis_->getDataNonConst(itmp).deepCopy(UnassembledSubspaceBases_[i]->getData(j)()); // Here, we copy data. Do we need to do this?
                        itmp++;
                    }
                }
            }
        }

        UnassembledBasesMaps_.resize(0);
        UnassembledSubspaceBases_.resize(0);

        UnassembledBasesMaps_.push_back(AssembledBasisMap_);
        UnassembledSubspaceBases_.push_back(AssembledBasis_);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::buildGlobalBasisMatrix(ConstXMapPtr rowMap,
                                                         ConstXMapPtr repeatedMap,
                                                         SC treshold)
    {
        FROSCH_ASSERT(!AssembledBasisMap_.is_null(),"AssembledBasisMap_.is_null().");
        FROSCH_ASSERT(!AssembledBasis_.is_null(),"AssembledBasis_.is_null().");

        GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,AssembledBasisMap_,AssembledBasisMap_->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!

        LO iD;
        SC valueTmp;
        GOVec indices;
        SCVec values;

        for (UN i=0; i<AssembledBasis_->getLocalLength(); i++) {
            indices.resize(0);
            values.resize(0);
            for (UN j=0; j<AssembledBasis_->getNumVectors(); j++) {
                valueTmp=AssembledBasis_->getData(j)[i];
                if (fabs(valueTmp)>treshold) {
                    indices.push_back( AssembledBasisMap_->getGlobalElement(j) );
                    values.push_back(valueTmp);
                }
            }
            iD = rowMap->getLocalElement(repeatedMap->getGlobalElement(i));

            if (iD!=-1) {
                GlobalBasisMatrix_->insertGlobalValues( repeatedMap->getGlobalElement(i) ,indices(),values());
            }
        }
        GlobalBasisMatrix_->fillComplete(AssembledBasisMap_,rowMap);
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::clearCoarseSpace()
    {
//        FROSCH_ASSERT(UnassembledBasesMaps_.size()>0,"UnassembledBasesMaps_.size()==0");
//        FROSCH_ASSERT(UnassembledSubspaceBases_.size()>0,"UnassembledSubspaceBases_.size()==0");

        UnassembledBasesMaps_.resize(0);
        UnassembledSubspaceBases_.resize(0);

        AssembledBasisMap_.reset();
        AssembledBasis_.reset();

        GlobalBasisMatrix_.reset();

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::checkForLinearDependencies()
    {
        FROSCH_ASSERT(false,"This is not implemented yet.");
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasBasisMap() const
    {
        return !AssembledBasisMap_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::XMapPtr CoarseSpace<SC,LO,GO,NO>::getBasisMap() const
    {
        FROSCH_ASSERT(!AssembledBasisMap_.is_null(),"AssembledBasisMap_.is_null().");
        return AssembledBasisMap_;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasAssembledBasis() const
    {
        return !AssembledBasis_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::XMultiVectorPtr CoarseSpace<SC,LO,GO,NO>::getAssembledBasis() const
    {
        FROSCH_ASSERT(!AssembledBasis_.is_null(),"AssembledBasis_.is_null().");
        return AssembledBasis_;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasGlobalBasisMatrix() const
    {
        return !GlobalBasisMatrix_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::XMatrixPtr CoarseSpace<SC,LO,GO,NO>::getGlobalBasisMatrix() const
    {
        FROSCH_ASSERT(!GlobalBasisMatrix_.is_null(),"GlobalBasisMatrix_.is_null().");
        return GlobalBasisMatrix_;
    }
}

#endif
