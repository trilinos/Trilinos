#ifndef MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP
#define MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP

#include "stk_search/BoxIdent.hpp"
#include "stk_search/HelperTraits.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Tree.hpp"
#include "stk_search/morton_lbvh/MortonLBVH_Search.hpp"
#include "stk_search/Box.hpp"
#include "stk_search/BoundingBox.hpp"
#include "stk_search/CommonSearchUtil.hpp"
#include "stk_search/DeviceMPIUtils.hpp"
#include "stk_search/HelperTraits.hpp"
#include "Kokkos_Core.hpp"
#include <vector>
#include <utility>

namespace stk::search {

template <typename DomainBoxType, typename DomainIdentProcType>
std::vector<Box<typename DomainBoxType::value_type>>
gather_all_processor_superset_domain_boxes(const std::vector<std::pair<DomainBoxType, DomainIdentProcType>> & localDomain,
                                           MPI_Comm & comm)
{
  std::vector<Box<typename DomainBoxType::value_type>> globalSupersetBoxes;

  Box<typename DomainBoxType::value_type> localSupersetBox;
  for (const auto & [box, ident] : localDomain) {
    stk::search::add_to_box(localSupersetBox, box);
  }

  stk::search::all_gather_helper(localSupersetBox, globalSupersetBoxes, comm);

  return globalSupersetBoxes;
}


namespace impl {
template <typename DomainView, typename ExecutionSpace>
class BoundingBoxReduction
{
  public:
    using DomainBoxType = typename DomainView::value_type::box_type;
    using Real          = typename DomainBoxType::value_type;
    using ResultBoxType = Box<Real>;

    using value_type = ResultBoxType;


    BoundingBoxReduction(DomainView localDomain) :
      m_localDomain(localDomain)
    {
      check_domain_or_range_view_parallel<DomainView, ExecutionSpace>();
    }

    ResultBoxType run(ExecutionSpace execSpace)
    {
      Kokkos::RangePolicy<ExecutionSpace> execPolicy(execSpace, size_t(0), m_localDomain.extent(0));

      ResultBoxType outputBox;
      Kokkos::parallel_reduce("local_bounding_box_reduction", execPolicy, *this, outputBox);
      execSpace.fence();

      return outputBox;
    }

    KOKKOS_INLINE_FUNCTION
    void init(ResultBoxType& val) const
    {
      constexpr Real min_val = Kokkos::Experimental::finite_min_v<Real>;
      constexpr Real max_val = Kokkos::Experimental::finite_max_v<Real>;

      val = ResultBoxType(max_val, max_val, max_val,
                          min_val, min_val, min_val);
    }

    KOKKOS_INLINE_FUNCTION
    void join(ResultBoxType& dest, const ResultBoxType& src) const
    {
      stk::search::add_to_box(dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int i, ResultBoxType& reductionBox) const
    {
      if constexpr (std::is_same_v<DomainBoxType, ResultBoxType>)
      {
        ResultBoxType box = m_localDomain(i).box;
        join(reductionBox, box);
      } else
      {
        DomainBoxType inputBox = m_localDomain(i).box;
        ResultBoxType outputBox(inputBox.get_x_min(), inputBox.get_y_min(), inputBox.get_z_min(),
                                inputBox.get_x_max(), inputBox.get_y_max(), inputBox.get_z_max());
        join(reductionBox, outputBox);
      }
    }

    private:
    DomainView m_localDomain;
};

}

template <typename DomainView, typename ExecutionSpace>
Kokkos::View<Box<typename DomainView::value_type::box_type::value_type>*, ExecutionSpace>
gather_all_processor_superset_domain_boxes(const DomainView & localDomain,
                                           ExecutionSpace execSpace,
                                           MPI_Comm comm)
{
  check_domain_or_range_view_parallel<DomainView, ExecutionSpace>();

  using BoundingBoxReductionType = impl::BoundingBoxReduction<DomainView, ExecutionSpace>;
  using BoxType = typename BoundingBoxReductionType::ResultBoxType;

  BoundingBoxReductionType func(localDomain);
  BoxType outputBox = func.run(execSpace);

  std::string name = "global_superset_boxes";
  Kokkos::View<BoxType*, ExecutionSpace> globalSupersetBoxes(Kokkos::view_alloc(name, Kokkos::WithoutInitializing), stk::parallel_machine_size(comm));
  stk::search::all_gather_helper(outputBox, globalSupersetBoxes, comm);

  return globalSupersetBoxes;
}


template<typename ExecutionSpace, typename DomainBoxType, typename DomainIdentProcType, typename RangeBoxType, typename RangeIdentProcType>
std::pair<std::vector<RangeBoxType>, std::vector<RangeIdentProcType>>
morton_extend_local_range_with_remote_boxes_that_might_intersect(
    const std::vector<std::pair<DomainBoxType, DomainIdentProcType>> & localDomain,
    const std::vector<std::pair<RangeBoxType, RangeIdentProcType>> & localRange,
    MPI_Comm comm,
    ExecutionSpace const& execSpace)
{
  using DomainValueType = typename DomainBoxType::value_type;
  using RangeValueType = typename RangeBoxType::value_type;

  const int numProcs = stk::parallel_machine_size(comm);
  const int procId = stk::parallel_machine_rank(comm);

  using StkDomainBoxType = stk::search::Box<DomainValueType>;
  using StkRangeBoxType = stk::search::Box<RangeValueType>;
  const std::vector<StkDomainBoxType> globalSupersetBoxes = gather_all_processor_superset_domain_boxes(localDomain, comm);

  using DomainViewType = Kokkos::View<BoxIdent<StkDomainBoxType,DomainIdentProcType>*,ExecutionSpace>;
  using RangeViewType = Kokkos::View<BoxIdent<StkRangeBoxType,RangeIdentProcType>*,ExecutionSpace>;
  using DomainTreeType = stk::search::MortonAabbTree<DomainViewType,ExecutionSpace>;
  using RangeTreeType = stk::search::MortonAabbTree<RangeViewType,ExecutionSpace>;

  DomainTreeType domainTree("Proc Domain Tree", localRange.size());
  RangeTreeType rangeTree("Proc Range Tree", globalSupersetBoxes.size());

  export_from_box_ident_proc_vec_to_morton_tree<DomainTreeType,RangeBoxType,RangeIdentProcType>(localRange, domainTree);
  export_from_box_vec_to_morton_tree<RangeTreeType,ExecutionSpace,StkDomainBoxType>(globalSupersetBoxes, rangeTree);
  domainTree.sync_to_device();
  rangeTree.sync_to_device();

  stk::search::CollisionList<ExecutionSpace> collisionList("Proc Collision List");
  stk::search::morton_lbvh_search<DomainTreeType, RangeTreeType, ExecutionSpace>(domainTree, rangeTree, collisionList, execSpace);
  collisionList.sync_from_device();

  using GlobalIdType = typename RangeIdentProcType::ident_type;
  using BoxIdPair = std::pair<RangeBoxType, GlobalIdType>;
  std::vector<std::vector<BoxIdPair>> sendList(numProcs);
  std::vector<std::vector<BoxIdPair>> recvList(numProcs);

  const unsigned numCollisions = collisionList.hm_idx();

  for (unsigned i = 0; i < numCollisions; ++i) {
    const int entityIndex = collisionList.hm_data(i, 0);
    const int remoteProcId = collisionList.hm_data(i, 1);
    const auto & [localRangeBox, localRangeIdentProc] = localRange[entityIndex];
    if (remoteProcId != procId) {
      sendList[remoteProcId].emplace_back(localRangeBox, localRangeIdentProc.id());
    }
  }

  stk::parallel_data_exchange_t(sendList, recvList, comm);

  std::pair<std::vector<RangeBoxType>, std::vector<RangeIdentProcType>> result;
  auto & [extendedRangeBoxes, remoteRangeIdentProcs] = result;

  extendedRangeBoxes.reserve(localRange.size());
  for (const auto & [box, identProc] : localRange) {
    extendedRangeBoxes.push_back(box);
  }

  for (size_t proc = 0; proc < recvList.size(); ++proc) {
    for (const auto & [box, id] : recvList[proc]) {
      extendedRangeBoxes.push_back(box);
      remoteRangeIdentProcs.emplace_back(id, proc);
    }
  }

  return result;
}

namespace impl {

template <typename RangeView, typename ExecutionSpace>
class FillGhostBoxBuffers
{
  public:
    using CollisionListType    = CollisionList<ExecutionSpace>;
    using SendDataType         = typename RangeView::value_type;
    using DeviceBufferAppender = impl::DeviceMPIBufferAppender<SendDataType, ExecutionSpace>;
    using DeviceBuffers        = typename DeviceBufferAppender::DeviceBuffers;

    FillGhostBoxBuffers(CollisionListType collisionList,
                        const RangeView& rangeBoxIdentProc,
                        ExecutionSpace execSpace,
                        MPI_Comm comm) :
      m_collisionList(collisionList),
      m_rangeBoxIdentProc(rangeBoxIdentProc),
      m_execSpace(execSpace),
      m_commRank(stk::parallel_machine_rank(comm)),
      m_deviceBufferAppender(stk::parallel_machine_size(comm), execSpace)
    {
      check_domain_or_range_view_parallel<RangeView, ExecutionSpace>();
    }

    DeviceBuffers run()
    {
      Kokkos::RangePolicy<ExecutionSpace> policy(m_execSpace, 0, m_collisionList.get_num_collisions());
      Kokkos::parallel_for("mpi_buffer_size_calc", policy, *this);
      m_execSpace.fence();

      m_deviceBufferAppender.allocate_buffers();

      Kokkos::parallel_for("mpi_buffer_fill", policy, *this);
      m_execSpace.fence();

      return m_deviceBufferAppender.getBuffers();
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(int idx) const
    {
      const int myProcId      = m_commRank;
      const int localBoxIdx   = m_collisionList.m_data(idx, 0);
      const int remoteProcId  = m_collisionList.m_data(idx, 1);

      if (remoteProcId != myProcId)
      {
        m_deviceBufferAppender.push_back(remoteProcId, m_rangeBoxIdentProc(localBoxIdx));
      }
    }

  private:
    CollisionListType m_collisionList;
    RangeView m_rangeBoxIdentProc;
    ExecutionSpace m_execSpace;
    int m_commRank;

    DeviceBufferAppender m_deviceBufferAppender;
};
}

template<typename DomainView, typename RangeView, typename ExecutionSpace>
std::pair<Kokkos::View<typename RangeView::value_type::box_type*, ExecutionSpace>,
          Kokkos::View<typename RangeView::value_type::ident_proc_type*, ExecutionSpace>>
morton_extend_local_range_with_remote_boxes_that_might_intersect(
    const DomainView & localDomain,
    const RangeView & localRange,
    ExecutionSpace execSpace,
    MPI_Comm comm)
{
  check_domain_or_range_view_parallel<DomainView, ExecutionSpace>();
  check_domain_or_range_view_parallel<RangeView, ExecutionSpace>();
  using DomainRealType               = typename DomainView::value_type::box_type::value_type;
  using RangeRealType               = typename RangeView::value_type::box_type::value_type;
  using RangeBoxType       = typename RangeView::value_type::box_type;
  using RangeIdentProcType = typename RangeView::value_type::ident_proc_type;
  using ViewType = Kokkos::View<Box<DomainRealType>*, ExecutionSpace>;
  using BoxIdentViewType = Kokkos::View<BoxIdent<Box<DomainRealType>,int>*,ExecutionSpace>;

  using MRangeBoxType = stk::search::Box<RangeRealType>;
  using MRangeViewType = Kokkos::View<BoxIdentProc<MRangeBoxType,RangeIdentProcType>*,ExecutionSpace>;
  using DomainTreeType = stk::search::MortonAabbTree<MRangeViewType, ExecutionSpace>;
  using RangeTreeType = stk::search::MortonAabbTree<BoxIdentViewType, ExecutionSpace>;

  using FillMPIBuffersType = impl::FillGhostBoxBuffers<RangeView, ExecutionSpace>;
  using DeviceBuffers = typename FillMPIBuffersType::DeviceBuffers;

  ViewType globalSupersetBoxes = gather_all_processor_superset_domain_boxes(localDomain, execSpace, comm);

  const bool setBoxesOnHost = false;
  DomainTreeType domainTree("Proc Domain Tree", localRange.size(), setBoxesOnHost);
  RangeTreeType rangeTree("Proc Range Tree", globalSupersetBoxes.size(), setBoxesOnHost);

  export_box_ident_view_to_morton_tree<RangeView,DomainTreeType,ExecutionSpace>(localRange, domainTree, execSpace);
  export_box_view_to_morton_tree<RangeTreeType,ExecutionSpace,ViewType>(globalSupersetBoxes, rangeTree, execSpace);
  execSpace.fence();
  domainTree.sync_to_device();
  rangeTree.sync_to_device();

  stk::search::CollisionList<ExecutionSpace> collisionList("Proc Collision List");
  stk::search::morton_lbvh_search<DomainTreeType, RangeTreeType, ExecutionSpace>(domainTree, rangeTree, collisionList, execSpace);

  FillMPIBuffersType fill_buffers(collisionList, localRange, execSpace, comm);
  DeviceBuffers deviceSendBuffers = fill_buffers.run();
  impl::DeviceDataExchangeUnknownPattern<typename DeviceBuffers::ValueType, ExecutionSpace> exchanger(deviceSendBuffers, execSpace, comm);
  DeviceBuffers deviceRecvBuffers = exchanger.communicate();

  size_t extendedRangeSize = localRange.size() + deviceRecvBuffers.buffers.extent(0);
  Kokkos::View<RangeBoxType*, ExecutionSpace> extendedRangeBoxes(Kokkos::ViewAllocateWithoutInitializing("extended_range_boxes"),
                                                                 extendedRangeSize);
  Kokkos::View<RangeIdentProcType*, ExecutionSpace> extendedRangeIdentProcs(Kokkos::ViewAllocateWithoutInitializing("extended_range_ident_procs"),
                                                                        deviceRecvBuffers.buffers.extent(0));

  Kokkos::RangePolicy<ExecutionSpace> policy(execSpace, 0, extendedRangeSize);
  size_t numLocalBoxes = localRange.extent(0);
  auto fillExtendedRange = KOKKOS_LAMBDA (size_t i)
  {
    if (i < numLocalBoxes)
    {
      extendedRangeBoxes(i) = localRange(i).box;
    } else
    {
      size_t bufferIdx = i - numLocalBoxes;
      extendedRangeBoxes(i) = deviceRecvBuffers.buffers(bufferIdx).box;
      extendedRangeIdentProcs(bufferIdx) = deviceRecvBuffers.buffers(bufferIdx).identProc;
    }
  };

  Kokkos::parallel_for("fill_extended_range", policy, fillExtendedRange);
  execSpace.fence();

  return std::make_pair(extendedRangeBoxes, extendedRangeIdentProcs);
}

}

#endif // MORTONLBVH_PARALLELCONSISTENCYUTILS_HPP
