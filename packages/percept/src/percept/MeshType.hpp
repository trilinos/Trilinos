// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef MeshType_hpp_
#define MeshType_hpp_

#include <percept/Percept.hpp>

#include <array>
#include <memory>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Bucket.hpp>

#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <percept/function/MDArray.hpp>
#include <percept/structured/StructuredCellIndex.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>

namespace percept {


enum NodeClassifyType {
    MS_VERTEX,
    MS_CURVE,
    MS_SURFACE,
    MS_VOLUME,
    MS_ON_BOUNDARY,
    MS_NOT_ON_BOUNDARY
};


#if defined(KOKKOS_ENABLE_CUDA)
            typedef double Double; //madbrew: long double's are problematic on gpu's as GPU's only use up to 64 bits for floating point representation. In some cases, they simply degrade into doubles but in others they cause illegal memory access errors for cuda
#else
            typedef long double Double;
//            typedef double Double;
#endif

	#ifdef KOKKOS_ENABLE_CUDA		
	  typedef Kokkos::CudaSpace     MemSpace ;
  	  typedef Kokkos::LayoutLeft    DataLayout;
  	  typedef Kokkos::Cuda          ExecSpace;
  	  //can be built either with serial or openmp options via the openmp=on option
#if defined(KOKKOS_ENABLE_OPENMP) // if you built with an nvidia compiler and have
                                // openmp=on in your build command
      typedef Kokkos::OpenMP     SecondaryExecSpace;
      typedef Kokkos::OpenMP     SecondaryMemSpace;
      typedef Kokkos::LayoutLeft SecondaryDataLayout;
#else
      typedef Kokkos::Serial        SecondaryExecSpace;
      typedef Kokkos::HostSpace     SecondaryMemSpace;
      typedef Kokkos::LayoutRight   SecondaryDataLayout;
#endif

#elif defined(KOKKOS_ENABLE_OPENMP)
          typedef Kokkos::OpenMP     ExecSpace;
  	  typedef Kokkos::OpenMP     MemSpace;
  	  typedef Kokkos::LayoutLeft DataLayout;

  	typedef Kokkos::OpenMP     SecondaryExecSpace;
  	typedef Kokkos::OpenMP     SecondaryMemSpace;
  	typedef Kokkos::LayoutLeft SecondaryDataLayout;

        #else
  	  typedef Kokkos::Serial      ExecSpace;
  	  typedef Kokkos::HostSpace   MemSpace;
  	  typedef Kokkos::LayoutRight DataLayout;

  	  typedef Kokkos::Serial        SecondaryExecSpace;
  	  typedef Kokkos::HostSpace     SecondaryMemSpace;
  	  typedef Kokkos::LayoutRight   SecondaryDataLayout;
	#endif
  	typedef Kokkos::View<double****, DataLayout , MemSpace > viewType;

    KOKKOS_INLINE_FUNCTION
    int device_safe_abs_int(int in)
    {
        return (in<0 ? (-1)*in : in );
    }

    KOKKOS_INLINE_FUNCTION
    int device_safe_sign(int value) { return value < 0 ? -1 : 1; }

    KOKKOS_INLINE_FUNCTION
    int device_safe_del(int v1, int v2) { return (int)(device_safe_abs_int(v1) == device_safe_abs_int(v2)); }

    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<int, 9> device_safe_transform_matrix(Kokkos::Array<int, 3>& trans_arr)
    {
        Kokkos::Array<int, 9> t_matrix;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          t_matrix[3 * i + j] = device_safe_sign(trans_arr[j]) * device_safe_del(trans_arr[j], i + 1);
        }
      }
      return t_matrix;
    }

    KOKKOS_INLINE_FUNCTION
    Kokkos::Array<int, 3> device_safe_transform_block_indices(Kokkos::Array<int, 3>& index_1,
                                    Kokkos::Array<int, 3>& trans_arr,
                                    Kokkos::Array<int, 3>& localRangeBeg,
                                    Kokkos::Array<int, 3>& donorRangeBeg)
    {
      Kokkos::Array<int, 9> t_matrix = device_safe_transform_matrix(trans_arr);

      Kokkos::Array<int, 3> diff;
      Kokkos::Array<int, 3> donor;

      diff[0] = index_1[0] - localRangeBeg[0];
      diff[1] = index_1[1] - localRangeBeg[1];
      diff[2] = index_1[2] - localRangeBeg[2];

      donor[0] =
          t_matrix[0] * diff[0] + t_matrix[1] * diff[1] + t_matrix[2] * diff[2] + donorRangeBeg[0];
      donor[1] =
          t_matrix[3] * diff[0] + t_matrix[4] * diff[1] + t_matrix[5] * diff[2] + donorRangeBeg[1];
      donor[2] =
          t_matrix[6] * diff[0] + t_matrix[7] * diff[1] + t_matrix[8] * diff[2] + donorRangeBeg[2];

      return donor;
    }

    template<typename T>
  	KOKKOS_INLINE_FUNCTION
  	T device_safe_max(T a, T b)
  	{
  	    return (a>b ? a : b);
  	}

    KOKKOS_INLINE_FUNCTION
    Double device_safe_abs_flt(Double in)
    {
        return (in<0 ? (-1)*in : in );
    }

//using StructuredCellIndex = std::array<unsigned,4>;  // i,j,k, block
class PerceptMesh;

struct SGridSelector {
    virtual bool operator()(StructuredCellIndex& indx) {
        VERIFY_MSG("not impl");
    }
};

struct SGridMeshGeometry {
    void normal_at(PerceptMesh *eMesh, StructuredCellIndex node,
            std::vector<double>& norm) {
        VERIFY_MSG("not impl");
    }

    int classify_node(StructuredCellIndex node_ptr,
            size_t& curveOrSurfaceEvaluator) {
        // FIXME
        return 2;
    }
};

struct MTSGridField {
    using Array4D = viewType;
    std::vector<std::shared_ptr<Array4D> > m_block_fields;
    std::string m_name;
    MTSGridField(const std::string& nm) :
            m_name(nm) {
    }
    //MTSGridField(const std::string& nm, StructuredBlock *sgrid) : m_name(mm) {}
    std::string name() {
        return m_name;
    }
};


struct SGridSizes {
  unsigned node_min[3];  // index min
  unsigned node_max[3];  // index max - all loops are for(i = node_min[0]; i <= node_max[0]; ++i)
  unsigned cell_min[3];
  unsigned cell_max[3];
  unsigned node_size[3];
  unsigned cell_size[3];

  unsigned node_size_global[3]; // this block may be parallel-distributed, this holds the size of the 'parent'/serial block
  //UInt pcell_size[3];
};

struct StructuredGrid {
    typedef SGridSelector MTSelector;
    typedef SGridMeshGeometry MTMeshGeometry;
    typedef StructuredCellIndex MTNode;
    typedef StructuredCellIndex MTElement;
    typedef StructuredCellIndex MTBucket;
#if STK_PERCEPT_LITE
    typedef double MTField;
#else
    typedef MTSGridField MTField;
#endif
    //    typedef int MTCellTopology;
    typedef CellTopologyData MTCellTopology;
    enum {
        NELEM_TYPES = 1
    };
};


 KOKKOS_INLINE_FUNCTION
 void sgrid_multi_dim_indices_from_index_node(const SGridSizes& block_sizes,
         const Kokkos::Array<unsigned int, 3> loop_orderings,
         const unsigned int& index,
         Kokkos::Array<unsigned, 3>& indx) {
     const int L0 = loop_orderings[0], L1 =
             loop_orderings[1], L2 = loop_orderings[2];
     const unsigned int sizes[3] = { 1 + block_sizes.node_max[L0]
             - block_sizes.node_min[L0], 1
             + block_sizes.node_max[L1]//SGridSizes block_sizes;
             - block_sizes.node_min[L1], 1//Kokkos::Array<unsigned int, 3> loop_orderings;
             + block_sizes.node_max[L2]
             - block_sizes.node_min[L2] };
     indx[L2] = block_sizes.node_min[L2]
             + (index / (sizes[L0] * sizes[L1]));
     indx[L1] = block_sizes.node_min[L1]
             + ((index / sizes[L0]) % sizes[L1]);
     indx[L0] = block_sizes.node_min[L0] + (index % sizes[L0]);

#if !defined(KOKKOS_ENABLE_CUDA) // exceptions cannot be called from the GPU
     unsigned int ii = indx[L0] - block_sizes.node_min[L0]
             + sizes[L0] * (indx[L1] - block_sizes.node_min[L1])
             + sizes[L0] * sizes[L1]
                     * (indx[L2] - block_sizes.node_min[L2]);
     VERIFY_OP_ON(ii, ==, index, "bad index");
#endif
 }


 KOKKOS_INLINE_FUNCTION
 void sgrid_multi_dim_indices_from_index_cell(const SGridSizes& block_sizes,
         const Kokkos::Array<unsigned int, 3> loop_orderings,
         const unsigned int& index,
         Kokkos::Array<unsigned, 3>& indx) {
     const int L0 = loop_orderings[0], L1 =
             loop_orderings[1], L2 = loop_orderings[2];
     const unsigned int sizes[3] = { 1 + block_sizes.cell_max[L0]
             - block_sizes.cell_min[L0], 1
             + block_sizes.cell_max[L1]//SGridSizes block_sizes;
             - block_sizes.cell_min[L1], 1//Kokkos::Array<unsigned int, 3> loop_orderings;
             + block_sizes.cell_max[L2]
             - block_sizes.cell_min[L2] };
     indx[L2] = block_sizes.cell_min[L2]
             + (index / (sizes[L0] * sizes[L1]));
     indx[L1] = block_sizes.cell_min[L1]
             + ((index / sizes[L0]) % sizes[L1]);
     indx[L0] = block_sizes.cell_min[L0] + (index % sizes[L0]);

#if !defined(KOKKOS_ENABLE_CUDA) // exceptions cannot be called from the GPU
     unsigned int ii = indx[L0] - block_sizes.cell_min[L0]
             + sizes[L0] * (indx[L1] - block_sizes.cell_min[L1])
             + sizes[L0] * sizes[L1]
                     * (indx[L2] - block_sizes.cell_min[L2]);
     VERIFY_OP_ON(ii, ==, index, "bad index");
#endif
 }


  class MeshGeometry;

  struct STKMesh {
    typedef stk::mesh::Selector MTSelector;
    typedef MeshGeometry MTMeshGeometry;
    typedef stk::mesh::Entity MTNode;
    typedef stk::mesh::Entity MTElement;
    typedef stk::mesh::Bucket MTBucket;
    typedef stk::mesh::FieldBase MTField;
    typedef CellTopologyData MTCellTopology;
    enum {NELEM_TYPES = 10 };
  };

  KOKKOS_INLINE_FUNCTION
     std::pair<bool,int> get_fixed_flag_sgrid(StructuredGrid::MTNode node_ptr, StructuredGrid::MTSelector *boundarySelector)
      {
     //   int dof = -1;
        std::pair<bool,int> ret(true,MS_VERTEX);
        //if the owner is something other than the top-level owner, the node
        // is on the boundary; otherwise, it isn't.
        bool& fixed = ret.first;
        int& type = ret.second;
        if (boundarySelector)
          {
            if ((*boundarySelector)(node_ptr))
              {
                fixed=true;
                type=MS_ON_BOUNDARY;
              }
            else
              {
                fixed=false;
                type = MS_NOT_ON_BOUNDARY;
              }
          }
        else
          {
     //#if defined(STK_PERCEPT_HAS_GEOMETRY)
     //       if (m_meshGeometry)
     //         {
     //           size_t curveOrSurfaceEvaluator;
     //           dof = m_meshGeometry->classify_node(node_ptr, curveOrSurfaceEvaluator);
     //           // vertex
     //           if (dof == 0)
     //             {
     //               fixed=true;
     //               type=MS_VERTEX;
     //             }
     //           // curve (for now we hold these fixed)
     //           else if (dof == 1)
     //             {
     //               fixed=true;
     //               type=MS_CURVE;
     //               //fixed=false;   // FIXME
     //             }
     //           // surface - also fixed
     //           else if (dof == 2)
     //             {
     //               //fixed=false;
     //               fixed=true;
     //               if (m_eMesh->get_smooth_surfaces())
     //                 {
     //                   fixed = false;
     //                 }
     //               type=MS_SURFACE;
     //               if (DEBUG_PRINT) std::cout << "tmp srk found surface node unfixed= " << node_ptr << std::endl;
     //             }
     //           // interior/volume - free to move
     //           else
     //             {
     //               fixed=false;
     //               type=MS_VOLUME;
     //             }
     //         }
     //       else
     //#endif
              {
                fixed=false;
                type=MS_VOLUME;
              }
          }
     //   if (DEBUG_PRINT) std::cout << "tmp srk classify node= " << node_ptr << " dof= " << dof << " fixed= " << fixed << " type= " << type << std::endl;

        return ret;
      }

  template<typename MeshType>
  unsigned get_num_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element);

  template<typename MeshType>
  const typename MeshType::MTNode *get_nodes(PerceptMesh *eMesh, typename MeshType::MTElement element, std::vector<typename MeshType::MTNode> *nodes );

  template<typename MeshType>
  bool MTisGhostNode(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  bool MTnode_locally_owned(PerceptMesh *m_eMesh, typename MeshType::MTNode node);

  template<typename MeshType>
  void MTcommFields(std::vector<const typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  template<typename MeshType>
  void MTsum_fields(std::vector<const typename MeshType::MTField*>& fields, PerceptMesh *m_eMesh);

  /// gets @param field data from @param node into @param fld
  template<typename MeshType>
  void get_field(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  template<typename MeshType>
  void get_field_new(double *fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// gets @param field data from @param node into @param fld[@param index]
  template<typename MeshType>
  void get_field(double *fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);

  /// sets @param field data from @param fld[@param index] into @param node
  template<typename MeshType>
  void set_field(const double * fld, unsigned size, int index, PerceptMesh *eMesh, typename MeshType::MTField *field, typename MeshType::MTNode node);


  template<typename T, size_t N>
  std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr)
  {
    for (unsigned i=0; i < N; ++i)
      os << arr[i] << (i == N-1 ? "" : " ");
    return os;
  }


}
#endif
