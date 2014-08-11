#ifndef SIERRA_SIERRA_MESH_DETAILS_BUCKET_LOCATION_HPP
#define SIERRA_SIERRA_MESH_DETAILS_BUCKET_LOCATION_HPP

#include <sierra/mesh/details/bucket_key.hpp>
#include <sierra/mesh/details/bucket_ordinal.hpp>

namespace sierra {
namespace mesh {
namespace details {

class bucket_location {
  public:
    bucket_location( const bucket_key & b = bucket_key(),
        const bucket_ordinal    & o = bucket_ordinal() )
      : m_bucket(b), m_ordinal(o)
    {}

    bool operator==(bucket_location rhs) const
    { return m_bucket == rhs.m_bucket  &&
      m_ordinal == rhs.m_ordinal;
    }

    bool operator!=(bucket_location rhs) const
    { return !(*this == rhs); }

    bool operator<(bucket_location rhs) const
    { return m_bucket < rhs.m_bucket  ||
      ( m_bucket == rhs.m_bucket && m_ordinal < rhs.m_ordinal);
    }

    bool operator<=(bucket_location rhs) const
    { return *this < rhs || *this == rhs;  }

    bool operator>(bucket_location rhs) const
    { return m_bucket > rhs.m_bucket  ||
      ( m_bucket == rhs.m_bucket && m_ordinal > rhs.m_ordinal);
    }

    bool operator>=(bucket_location rhs) const
    { return *this > rhs || *this == rhs;  }

    const bucket_key     &  bucket()  const { return m_bucket; }
    const bucket_ordinal &  ordinal() const { return m_ordinal; }

  private:
    bucket_key        m_bucket;
    bucket_ordinal    m_ordinal;
};

inline bucket_key bucket(bucket_location e) { return e.bucket(); }
inline bucket_ordinal ordinal(bucket_location e) { return e.ordinal(); }

} // details
} // mesh
} // sierra

#endif //SIERRA_SIERRA_MESH_DETAILS_BUCKET_LOCATION_HPP
