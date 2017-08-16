#ifndef _STKMESHADAPTERFORZOLTAN2_HPP_
#define _STKMESHADAPTERFORZOLTAN2_HPP_

#include <stddef.h>                     // for size_t, NULL
#include <sys/types.h>                  // for ssize_t
#include <Zoltan2_MeshAdapter.hpp>      // for MeshEntityType, etc
#include <ostream>                      // for basic_ostream, char_traits, etc
#include <vector>                       // for vector
#include "Zoltan2_Adapter.hpp"
#include "Zoltan2_InputTraits.hpp"      // for BasicUserTypes
#include "balanceTypes.hpp"
#include "Zoltan2ParallelGraph.hpp"

typedef Zoltan2::BasicUserTypes<double, BalanceLocalNumber, BalanceGlobalNumber> stkdata_t;

class StkMeshZoltanAdapter : public Zoltan2::MeshAdapter<stkdata_t>
{
public:
    size_t getGlobalNumOf(Zoltan2::MeshEntityType etype) const;

    typedef Zoltan2::MeshAdapter<stkdata_t> base_adapter_t;

    StkMeshZoltanAdapter(const Zoltan2ParallelGraph &graph);

    virtual ~StkMeshZoltanAdapter() { }

    virtual size_t getLocalNumOf(Zoltan2::MeshEntityType etype) const;
    virtual void getIDsViewOf(Zoltan2::MeshEntityType etype, BalanceGlobalNumber const *&Ids) const;

    virtual int getDimension() const;

    virtual void getCoordinatesViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&coords, int &stride, int coordDim) const;

    virtual int getNumWeightsPerOf(Zoltan2::MeshEntityType etype) const;

    virtual void getWeightsViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&weights, int &stride, int idx = 0) const;

    virtual bool avail2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const;

    virtual size_t getLocalNum2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const;

    virtual void get2ndAdjsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *&adjacencyIds) const;

    virtual int getNumWeightsPer2ndAdj(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const;

    virtual void get2ndAdjWeightsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const scalar_t *&weights, int &stride, int idx) const;

    void debuggingInfo(int proc_id,  std::ofstream& out) const;

private:
    const Zoltan2ParallelGraph &mGraph;

public: // defaultish
    virtual bool availAdjs(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target) const
    {
        return false;
    }

    virtual size_t getLocalNumAdjs(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target) const
    {
        return 0;
    }

    virtual void getAdjsView(Zoltan2::MeshEntityType source, Zoltan2::MeshEntityType target, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *& adjacencyIds) const
    {
        offsets = NULL;
        adjacencyIds = NULL;
        Z2_THROW_NOT_IMPLEMENTED
    }

    virtual bool useDegreeAsWeightOf(Zoltan2::MeshEntityType etype, int idx) const
    {
        return false;
    }
};

#endif
