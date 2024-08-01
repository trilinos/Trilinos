// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ioss_Utils.h>
#include <catalyst_tests/Iocatalyst_BlockMesh.h>

namespace Iocatalyst {

  BlockMesh::ID BlockMesh::_id = 1;

  BlockMesh::BlockMesh()
  {
    partition.id       = 0;
    partition.size     = 1;
    extents.i          = 1;
    extents.j          = 1;
    extents.k          = 1;
    partitionExtents.i = 1;
    partitionExtents.j = 1;
    partitionExtents.k = 1;
    partitionStart.i   = 0;
    partitionStart.j   = 0;
    partitionStart.k   = 0;
    id                 = _id++;
  }

  BlockMesh::~BlockMesh() {}

  void BlockMesh::init(const Partition &part, const Extent &numBlocks, const Extent &origin)
  {
    if (part.id < 0 || part.size < 0 || part.id >= part.size) {
      std::ostringstream errmsg;
      errmsg << "Invalid partition: id = " << part.id << std::string(", size = ") << part.size
             << "\n";
      IOSS_ERROR(errmsg);
    }
    this->partition = part;
    this->extents   = numBlocks;
    this->origin    = origin;
    splitBlock();
  }

  void BlockMesh::splitBlock()
  {
    // Split algorithm from vtkExtentTranslator.cxx SplitExtent()

    unsigned long size[3];
    int           numPiecesInFirstHalf;
    int           splitAxis;
    long int      mid;
    int           ext[6];
    int           numPieces = getPartition().size;
    int           piece     = getPartition().id;
    fillExtents(ext);
    setPartitionFromExtents(ext);

    while (numPieces > 1) {
      size[0] = ext[1] - ext[0];
      size[1] = ext[3] - ext[2];
      size[2] = ext[5] - ext[4];

      if (size[2] >= size[1] && size[2] >= size[0] && size[2] / 2 >= 1) {
        splitAxis = 2;
      }
      else if (size[1] >= size[0] && size[1] / 2 >= 1) {
        splitAxis = 1;
      }
      else if (size[0] / 2 >= 1) {
        splitAxis = 0;
      }
      else {
        splitAxis = -1;
      }

      if (splitAxis == -1) {
        if (piece == 0) {
          numPieces = 1;
        }
        else {
          setPartitionEmpty();
          return;
        }
      }
      else {
        numPiecesInFirstHalf = (numPieces / 2);
        mid                  = size[splitAxis];
        mid                  = (mid * numPiecesInFirstHalf) / numPieces + ext[splitAxis * 2];
        if (piece < numPiecesInFirstHalf) {
          ext[splitAxis * 2 + 1] = mid;
          numPieces              = numPiecesInFirstHalf;
        }
        else {
          ext[splitAxis * 2] = mid;
          numPieces          = numPieces - numPiecesInFirstHalf;
          piece -= numPiecesInFirstHalf;
        }
      }
    }
    setPartitionFromExtents(ext);
  }

  void BlockMesh::fillExtents(int *ext)
  {
    if (getExtents().i == 0) {
      ext[0] = 0;
      ext[1] = -1;
    }
    else {
      ext[0] = 0;
      ext[1] = getExtents().i;
    }

    if (getExtents().j == 0) {
      ext[2] = 0;
      ext[3] = -1;
    }
    else {
      ext[2] = 0;
      ext[3] = getExtents().j;
    }

    if (getExtents().k == 0) {
      ext[4] = 0;
      ext[5] = -1;
    }
    else {
      ext[4] = 0;
      ext[5] = getExtents().k;
    }
  }

  void BlockMesh::setPartitionFromExtents(int ext[6])
  {
    int sizeX = ext[1] - ext[0];
    int sizeY = ext[3] - ext[2];
    int sizeZ = ext[5] - ext[4];
    if (sizeX <= 0 || sizeY <= 0 || sizeZ <= 0) {
      setPartitionEmpty();
      return;
    }

    partitionExtents.i = sizeX;
    partitionStart.i   = ext[0];

    partitionExtents.j = sizeY;
    partitionStart.j   = ext[2];

    partitionExtents.k = sizeZ;
    partitionStart.k   = ext[4];
  }

  void BlockMesh::setPartitionEmpty()
  {
    partitionExtents.i = 0;
    partitionExtents.j = 0;
    partitionExtents.k = 0;
    partitionStart.i   = 0;
    partitionStart.j   = 0;
    partitionStart.k   = 0;
  }

  bool BlockMesh::isPartitionEmpty() const
  {
    return partitionExtents.i == 0 || partitionExtents.j == 0 || partitionExtents.k == 0;
  }

  BlockMesh::ID BlockMesh::getID() const { return id; }

  int BlockMesh::getNumPoints() const { return getNumInBlockMesh(POINT_OFFSET); }

  int BlockMesh::getNumPartitionPoints() const { return getNumInPartition(POINT_OFFSET); }

  BlockMesh::IDList BlockMesh::getPartitionPointIDs() const
  {
    return getPartitionIDs(POINT_OFFSET);
  }

  BlockMesh::ID BlockMesh::getPointIDfromCoords(unsigned int i, unsigned int j,
                                                unsigned int k) const
  {
    Extent coords = {i, j, k};
    Extent bounds = {extents.i + POINT_OFFSET, extents.j + POINT_OFFSET, extents.k + POINT_OFFSET};
    return getIDfromCoords(coords, bounds);
  }

  BlockMesh::Point BlockMesh::getPointCoordsForPointID(ID pointID) const
  {
    Point  p;
    Extent bounds = {extents.i + POINT_OFFSET, extents.j + POINT_OFFSET, extents.k + POINT_OFFSET};
    Extent ext    = getCoordsForID(pointID, bounds);
    p.x           = origin.i + ext.i * BLOCK_LENGTH;
    p.y           = origin.j + ext.j * BLOCK_LENGTH;
    p.z           = origin.k + ext.k * BLOCK_LENGTH;
    return p;
  }

  BlockMesh::ID BlockMesh::getGlobalIDForPointID(ID pointID)
  {
    Extent pointExt = getExtents();
    pointExt.i += 1;
    pointExt.j += 1;
    pointExt.k += 1;
    auto coords = getCoordsForID(pointID, pointExt);
    coords.i += origin.i;
    coords.j += origin.j;
    coords.k += origin.k;
    return getIDfromCoords(coords, getGlobalPointExtents());
  }

  int BlockMesh::getNumBlocks() const { return getNumInBlockMesh(BLOCK_OFFSET); }

  int BlockMesh::getNumPartitionBlocks() const { return getNumInPartition(BLOCK_OFFSET); }

  BlockMesh::IDList BlockMesh::getPartitionBlockIDs() const
  {
    return getPartitionIDs(BLOCK_OFFSET);
  }

  BlockMesh::BlockConn BlockMesh::getBlockConnectivityPointIDs(ID blockID) const
  {
    BlockConn conn;
    Extent bounds = {extents.i + BLOCK_OFFSET, extents.j + BLOCK_OFFSET, extents.k + BLOCK_OFFSET};
    Extent e      = getCoordsForID(blockID, bounds);
    conn[0]       = getPointIDfromCoords(e.i, e.j, e.k);
    conn[1]       = getPointIDfromCoords(e.i + 1, e.j, e.k);
    conn[2]       = getPointIDfromCoords(e.i + 1, e.j + 1, e.k);
    conn[3]       = getPointIDfromCoords(e.i, e.j + 1, e.k);
    conn[4]       = getPointIDfromCoords(e.i, e.j, e.k + 1);
    conn[5]       = getPointIDfromCoords(e.i + 1, e.j, e.k + 1);
    conn[6]       = getPointIDfromCoords(e.i + 1, e.j + 1, e.k + 1);
    conn[7]       = getPointIDfromCoords(e.i, e.j + 1, e.k + 1);
    return conn;
  }

  BlockMesh::ID BlockMesh::getGlobalIDForBlockID(ID blockID)
  {
    auto coords = getCoordsForID(blockID, getExtents());
    coords.i += origin.i;
    coords.j += origin.j;
    coords.k += origin.k;
    return getIDfromCoords(coords, getGlobalBlockExtents());
  }

  BlockMesh::IDList BlockMesh::getPartitionIDs(unsigned int offset) const
  {
    BlockMesh::IDList ids;
    Extent            bounds = {extents.i + offset, extents.j + offset, extents.k + offset};
    if (!isPartitionEmpty()) {
      for (unsigned int k = partitionStart.k; k < partitionStart.k + partitionExtents.k + offset;
           k++) {
        for (unsigned int j = partitionStart.j; j < partitionStart.j + partitionExtents.j + offset;
             j++) {
          for (unsigned int i = partitionStart.i;
               i < partitionStart.i + partitionExtents.i + offset; i++) {
            Extent coords = {i, j, k};
            ids.push_back(getIDfromCoords(coords, bounds));
          }
        }
      }
    }
    return ids;
  }

  BlockMesh::ID BlockMesh::getIDfromCoords(Extent coords, Extent bounds)
  {
    return coords.k * bounds.i * bounds.j + coords.j * bounds.i + coords.i + 1;
  }

  BlockMesh::Extent BlockMesh::getCoordsForID(ID id, Extent bounds)
  {
    Extent       ext;
    int          zeroBasedID = id - 1;
    unsigned int sizeXY      = bounds.i * bounds.j;
    ext.k                    = zeroBasedID / sizeXY;
    ext.j                    = (zeroBasedID - ext.k * sizeXY) / bounds.i;
    ext.i                    = zeroBasedID - ext.k * sizeXY - ext.j * bounds.i;
    return ext;
  }

  int BlockMesh::getNumInPartition(unsigned int offset) const
  {
    return (partitionExtents.i + offset) * (partitionExtents.j + offset) *
           (partitionExtents.k + offset);
  }

  int BlockMesh::getNumInBlockMesh(unsigned int offset) const
  {
    return (extents.i + offset) * (extents.j + offset) * (extents.k + offset);
  }

  std::map<std::string, double>* BlockMesh::getTransientCellFieldMap()
  {
    return &(this->transientCellFields);
  }

  std::map<std::string, double>* BlockMesh::getTransientPointFieldMap()
  {
    return &(this->transientPointFields);
  }

  void BlockMesh::addTransientCellField(std::string f_name, double f_value)
  {
    this->transientCellFields.insert({ f_name, f_value });
  }

  void BlockMesh::addTransientPointField(std::string f_name, double f_value)
  {
    this->transientPointFields.insert({ f_name, f_value });
  }

} // namespace Iocatalyst
