#ifndef stk_util_parallel_MPITag
#define stk_util_parallel_MPITag

#include <memory>
#include <ostream>
#include "stk_util/parallel/Parallel.hpp"   // for MPI
#include "stk_util/parallel/MPICommKey.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire


namespace stk {

class MPITagManager;

namespace impl {

class MPITagData
{
  public:

    MPITagData(MPITagManager* manager,
               MPI_Comm comm,
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
               MPI_Comm commInternal,
#endif
               int tag) :
      m_manager(manager),
      m_comm(comm),
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
      m_commInternal(commInternal),
#endif
      m_tag(tag)
    {}

    ~MPITagData();

    MPITagData(const MPITagData&) = delete;

    MPITagData& operator=(const MPITagData&) = delete;

    int get_tag() { return m_tag;}

    MPI_Comm get_comm() { return m_comm;}

#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    MPI_Comm get_comm_internal() { return m_commInternal; }
#endif 

    MPIKeyManager::CommKey get_comm_key();

    void set_free() { m_isFree = true;}

    bool is_free() { return m_isFree; }

    MPITagManager* get_tag_manager() { return m_manager; }

  private:
    MPITagManager* m_manager;
    MPI_Comm m_comm;
#ifdef MPI_KEY_MANAGER_COMM_DESTRUCTOR_CALLBACK_BROKEN
    MPI_Comm m_commInternal;
#endif
    int m_tag;
    bool m_isFree = false;
};

}  // namespace


// an MPI tag that is currently in use.  Note that it is implicitly convertable
// to int, so it can be passed directly into MPI routines.
class MPITag
{
  public:
    MPITag() = default;

    explicit MPITag(std::shared_ptr<impl::MPITagData> data) :
      m_data(data)
    {}

    operator int() const {
      STK_ThrowRequireMsg(m_data, "Cannot convert unitialized MPITag to int");
      return m_data->get_tag();
    }

    MPI_Comm get_comm() const {
      STK_ThrowRequireMsg(m_data, "Cannot call get_comm on unitialized MPITag");
      return m_data->get_comm();
    }

  private:
    void set_free()
    {
      if(m_data) {
        m_data->set_free();
      }
    }

    std::shared_ptr<impl::MPITagData> m_data;

    friend MPITagManager;

    friend bool operator==(const MPITag& lhs, const MPITag& rhs);
    friend bool operator!=(const MPITag& lhs, const MPITag& rhs);
    friend std::ostream& operator<<(std::ostream& os, const MPITag& tag);
};

inline bool operator==(const MPITag& lhs, const MPITag& rhs)
{
  STK_ThrowRequireMsg(lhs.m_data->get_tag_manager() == rhs.m_data->get_tag_manager(),
                  "Cannot compare MPITags on different MPITagManagers");

  return static_cast<int>(lhs) == static_cast<int>(rhs) &&
         !lhs.m_data->is_free() && !rhs.m_data->is_free() &&
         lhs.m_data->get_comm_key() == rhs.m_data->get_comm_key();
}

inline bool operator!=(const MPITag& lhs, const MPITag& rhs)
{
  return !(lhs == rhs);
}

inline std::ostream&  operator<<(std::ostream& os, const MPITag& tag)
{
  os << "MPITag with value " << tag.m_data->get_tag() << " on Comm with key "
     << tag.m_data->get_comm_key() << ", is_free = " << tag.m_data->is_free();

  return os;
}

} // namespace

#endif
