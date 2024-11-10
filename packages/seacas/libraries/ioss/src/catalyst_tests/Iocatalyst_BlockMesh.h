// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#pragma once

#include "iocatalyst_export.h"
#include <array>
#include <vector>

namespace Iocatalyst {

  class IOCATALYST_EXPORT BlockMesh
  {
  public:
    struct Partition
    {
      int id;
      int size;
    };

    struct Point
    {
      double x;
      double y;
      double z;
    };

    struct Extent
    {
      unsigned int i;
      unsigned int j;
      unsigned int k;
    };

    using BlockConn                        = std::array<int, 8>;
    using IDList                           = std::vector<int>;
    using ID                               = unsigned int;
    static const unsigned int BLOCK_OFFSET = 0;
    static const unsigned int POINT_OFFSET = 1;
    static constexpr double   BLOCK_LENGTH = 1.0;

    static const unsigned int I_GLOBAL = 1000;
    static const unsigned int J_GLOBAL = 1000;
    static const unsigned int K_GLOBAL = 1000;

    BlockMesh();
    ~BlockMesh();

    void init(const Partition &part, const Extent &numBlocks, const Extent &origin);

    const Partition &getPartition() const { return partition; }

    const Extent &getOrigin() const { return origin; }

    const Extent &getExtents() const { return extents; }

    const Extent &getPartitionExtents() const { return partitionExtents; }
    const Extent &getPartitionStart() const { return partitionStart; }

    ID getID() const;

    bool isPartitionEmpty() const;

    int           getNumBlocks() const;
    int           getNumPartitionBlocks() const;
    IDList        getPartitionBlockIDs() const;
    BlockConn     getBlockConnectivityPointIDs(ID blockID) const;
    static Extent getGlobalBlockExtents() { return {I_GLOBAL, J_GLOBAL, K_GLOBAL}; };
    ID            getGlobalIDForBlockID(ID blockID);

    int           getNumPoints() const;
    int           getNumPartitionPoints() const;
    IDList        getPartitionPointIDs() const;
    ID            getPointIDfromCoords(unsigned int i, unsigned int j, unsigned int k) const;
    Point         getPointCoordsForPointID(ID pointID) const;
    static Extent getGlobalPointExtents() { return {I_GLOBAL + 1, J_GLOBAL + 1, K_GLOBAL + 1}; };
    ID            getGlobalIDForPointID(ID pointID);

    static Extent getCoordsForID(ID id, Extent bounds);
    static ID     getIDfromCoords(Extent coords, Extent bounds);

    void          addTransientCellField(std::string f_name, double f_value);
    void          addTransientPointField(std::string f_name, double f_value);

    std::map<std::string, double>* getTransientCellFieldMap();
    std::map<std::string, double>* getTransientPointFieldMap();

  private:
    Partition partition;
    Extent    origin;
    Extent    extents;
    Extent    partitionExtents;
    Extent    partitionStart;
    Extent    globalBlockExtents;
    void      splitBlock();
    void      fillExtents(int *ext);
    void      setPartitionFromExtents(int ext[6]);
    void      setPartitionEmpty();
    IDList    getPartitionIDs(unsigned int offset) const;
    int       getNumInPartition(unsigned int offset) const;
    int       getNumInBlockMesh(unsigned int offset) const;
    ID        id;
    static ID _id;

    std::map<std::string, double> transientCellFields;
    std::map<std::string, double> transientPointFields;
  };

} // namespace Iocatalyst
