#ifndef stk_util_parallel_IAllreduceReplacement
#define stk_util_parallel_IAllreduceReplacement

#include <stddef.h>
#include <functional>
#include "Parallel.hpp"
#include "TreeReductionIndexer.hpp"

namespace stk {
namespace impl {

//TODO: Move the ifdef in USE_IBARRIER_REPLACEMENT in DeletionGroup.cpp
//      to here, and make this class call MPI_Iallreduce whenever possible
//      or maybe just do this for IbarrierReplacement?

template <typename T>
class IAllreduceReplacement
{
  enum class ReductionState
  {
    NotStarted,
    ReductionUp,
    FinishReductionUp,
    ReductionDown,
    FinishReductionDown,
    Completed
  };

  public:
    using Op = std::function<T(const T&, const T&)>;

    explicit IAllreduceReplacement(MPI_Comm comm, MPI_Datatype datatype, int tag, int reductionFactor=2) :
      m_comm(comm),
      m_tag(tag),
      m_reductionFactor(reductionFactor),
      m_indexer(stk::parallel_machine_size(comm), reductionFactor),
      m_myrank(stk::parallel_machine_rank(comm)),
      m_childRequests(reductionFactor-1, MPI_REQUEST_NULL),
      m_childRanks(reductionFactor-1),
      m_datatype(datatype),
      m_recv_buffers(reductionFactor-1)
    {}

    void startReduction(Op op, const std::vector<T>& input_data, std::vector<T>& output_data)
    {
      STK_ThrowRequireMsg(m_state == ReductionState::NotStarted || m_state == ReductionState::Completed,
                      "Cannot start a new barrier while the old one is still in progress");

      setup(op, input_data, output_data);

      if (m_indexer.get_tree_depth() == 1)
      {
        m_state = ReductionState::Completed;
        return;
      }

      bool isLastStep = start_reduction_up();
      m_state = isLastStep ? ReductionState::FinishReductionUp : ReductionState::ReductionUp;
    }

    bool progressReduction()
    {
      if (m_state == ReductionState::Completed)
      {
        return true;
      }

      switch (m_state)
      {
        case ReductionState::ReductionUp:
        {
          bool isFinished = reduction_up_step();
          if (isFinished)
          {
            m_state = ReductionState::FinishReductionUp;
          }
          break;
        }

        case ReductionState::FinishReductionUp:
        {
          bool isFinished = finish_reduction_up();
          if (isFinished)
          {
            start_reduction_down();
            m_state = ReductionState::ReductionDown;
          }

          break;
        }

        case ReductionState::ReductionDown:
        {
          bool isFinished = reduction_down_step();
          if (isFinished)
          {
            m_state = ReductionState::FinishReductionDown;
          }

          break;
        }

        case ReductionState::FinishReductionDown:
        {
          bool isFinished = finish_reduction_down();
          if (isFinished)
          {
            m_state = ReductionState::Completed;
          }
          break;
        }

        default:
          STK_ThrowRequireMsg(false, "Reduction not started, cannot progress it");
      };

      return m_state == ReductionState::Completed;
    }

    void finishReduction()
    {
      while(m_state != ReductionState::Completed)
      {
        progressReduction();
      }
    }

  private:

    void setup(Op op, const std::vector<T>& input_data, std::vector<T>& output_data)
    {
      STK_ThrowRequireMsg(input_data.size() == output_data.size(), "input and output buffers must be same size");

      m_op = op;
      m_output_buffer = &(output_data);

      for (int i=0; i < m_reductionFactor-1; ++i)
      {
        m_recv_buffers[i].resize(input_data.size());
      }

      for (size_t i=0; i < input_data.size(); ++i)
      {
        output_data[i] = input_data[i];
      }
    }

    //-------------------------------------------------------------------------
    // reduction up

    bool start_reduction_up()
    {
      bool isMyLastStep = reduction_up_step_impl();

      return isMyLastStep;
    }

    bool reduction_up_step()
    {
      bool isMyLastStep = false;
      int flag;
      MPI_Testall(m_childRequests.size(), m_childRequests.data(), &flag, MPI_STATUSES_IGNORE);
      if (flag)
      {
        apply_reduction();
        isMyLastStep = reduction_up_step_impl();
      }


      return isMyLastStep;
    }

    bool finish_reduction_up()
    {
      int isFinished = true;
      if (m_indexer.is_rank_on_level(m_indexer.get_tree_depth() - 1, m_myrank))
      {
        MPI_Testall(m_childRequests.size(), m_childRequests.data(), &isFinished, MPI_STATUSES_IGNORE);
        if (isFinished)
        {
          apply_reduction();
        }
      }

      return isFinished;
    }

    bool reduction_up_step_impl()
    {
      bool isMyLastStep = m_currentLevel + 1 == m_indexer.get_tree_depth() - 1 ||
                          !(m_indexer.is_rank_on_level(m_currentLevel + 1, m_myrank));

      if (m_indexer.is_rank_on_level(m_currentLevel + 1, m_myrank))
      {        
        m_currentLevel++;
        start_recv_from_children();
      }  else
      {
        start_send_to_parent();
      }

      return isMyLastStep;
    }

    void apply_reduction()
    {
      auto& output_buffer = *m_output_buffer;
      for (size_t j=0; j < output_buffer.size(); ++j)
      {
        for (int i=0; i < m_numChildren; ++i)
        {
          output_buffer[j] = m_op(output_buffer[j], m_recv_buffers[i][j]);
        }
      }
    }


    //-------------------------------------------------------------------------
    // reduction down

    void start_reduction_down()
    {
      if (m_indexer.is_rank_on_level(m_indexer.get_tree_depth() - 1, m_myrank))
      {
          start_send_to_children();
          m_currentLevel--;
      } else
      {
        start_recv_from_parent();
      }
    }

    bool reduction_down_step()
    {
      if (m_currentLevel > 0)
      {
        int flag;
        MPI_Testall(m_childRequests.size(), m_childRequests.data(), &flag, MPI_STATUSES_IGNORE);
        if (flag)
        {
          MPI_Test(&m_recvFromParentRequest, &flag, MPI_STATUS_IGNORE);
          if (flag)
          {
            start_send_to_children();
            m_currentLevel--;
          }
        }
      }

      bool isFinished = m_currentLevel == 0;
      return isFinished;
    }

    bool finish_reduction_down()
    {
      int flag1, flag2, flag3;
      MPI_Testall(m_childRequests.size(), m_childRequests.data(), &flag1, MPI_STATUSES_IGNORE); 

      MPI_Test(&m_recvFromParentRequest, &flag2, MPI_STATUS_IGNORE);
      MPI_Test(&m_sendToParentRequest, &flag3, MPI_STATUS_IGNORE);

      return flag1 && flag2 && flag3;  
    }


    //-------------------------------------------------------------------------
    // helper functions
    void start_recv_from_children()
    {
      m_numChildren = m_indexer.get_child_ranks(m_currentLevel-1, m_myrank, m_childRanks);
      for (int i=0; i < m_numChildren; ++i)
      {
        MPI_Irecv(m_recv_buffers[i].data(), m_recv_buffers[i].size(), m_datatype,
                  m_childRanks[i], m_tag, m_comm, &(m_childRequests[i]));
      }      
    }

    void start_send_to_parent()
    {
      int parentRank = m_indexer.get_parent_rank(m_currentLevel+1, m_myrank);
      MPI_Isend(m_output_buffer->data(), m_output_buffer->size(), m_datatype, parentRank, 
                m_tag, m_comm, &(m_sendToParentRequest));    
    }

    void start_send_to_children()
    {
      int numChildren = m_indexer.get_child_ranks(m_currentLevel-1, m_myrank, m_childRanks);
      for (int i=0; i < numChildren; ++i)
      {
        MPI_Isend(m_output_buffer->data(), m_output_buffer->size(), m_datatype, 
                  m_childRanks[i], m_tag, m_comm, &(m_childRequests[i]));
      }
    }

    void start_recv_from_parent()
    {
      int parentRank = m_indexer.get_parent_rank(m_currentLevel+1, m_myrank);
      MPI_Irecv(m_output_buffer->data(), m_output_buffer->size(), m_datatype, 
                parentRank, m_tag, m_comm, &(m_recvFromParentRequest));
    }

    MPI_Comm m_comm;
    int m_tag;
    int m_reductionFactor;
    TreeReductionIndexer m_indexer;
    int m_myrank;
    int m_currentLevel = 0;
    int m_numChildren = 0;
    std::vector<MPI_Request> m_childRequests;
    std::vector<int> m_childRanks;
    MPI_Request m_sendToParentRequest   = MPI_REQUEST_NULL;
    MPI_Request m_recvFromParentRequest = MPI_REQUEST_NULL;

    Op m_op;
    MPI_Datatype m_datatype;
    std::vector<T>* m_output_buffer;
    std::vector< std::vector<T> > m_recv_buffers;
    ReductionState m_state = ReductionState::NotStarted;
};

}
}

#endif
