#ifndef stk_search_BoundingBoxCompare_hpp
#define stk_search_BoundingBoxCompare_hpp

namespace stk {
namespace search {
namespace box {
namespace compare{

enum {
LOWER = 0,
MIDDLE = 1,
UPPER = 2,
LENGTH = 3,
KEY = 4
};

template <class BoundingBox, int Mode>
class Compare;

template <class BoundingBox, int Mode>
class Partition;

template<class BoundingBox>
class Compare<BoundingBox, LOWER> {
  public:
    Compare(int axis) : m_axis(axis) {}
    inline bool operator ()(const BoundingBox &a, const BoundingBox &b) const {
      return a.lower(m_axis) < b.lower(m_axis);
    }
  private:
    const int m_axis;
};

template<class BoundingBox>
class Compare<BoundingBox, UPPER> {
  public:
    Compare(int axis) : m_axis(axis) {}
    inline bool operator ()(const BoundingBox &a, const BoundingBox &b) const {
      return a.upper(m_axis) < b.upper(m_axis);
    }
  private:
    const int m_axis;
};

template<class BoundingBox>
class Compare<BoundingBox, MIDDLE> {
  public:
    Compare(int axis) : m_axis(axis) {}
    inline bool operator ()(const BoundingBox &a, const BoundingBox &b) const {
      return a.middle(m_axis) < b.middle(m_axis);
    }
  private:
    const int m_axis;
};

template<class BoundingBox>
class Compare<BoundingBox, LENGTH> {
  public:
    Compare(int axis) : m_axis(axis) {}
    inline bool operator ()(const BoundingBox &a, const BoundingBox &b) const {
      return a.length(m_axis) < b.length(m_axis);
    }
  private:
    const int m_axis;
};

template<class BoundingBox>
class Compare<BoundingBox, KEY> {
  public:
    Compare(int /*axis*/ = 0) {}
    inline bool operator ()(const BoundingBox &a, const BoundingBox &b) const {
      return a.key < b.key;
    }
  private:
};

template<class BoundingBox>
class Partition<BoundingBox, LOWER> {
  public:
    Partition(int axis, typename BoundingBox::Data split_plane) : m_axis(axis), m_split_plane(split_plane) {}
    inline bool operator ()(const BoundingBox &a) const {
      return a.lower(m_axis) < m_split_plane;
    }
  private:
    const int m_axis;
    typename BoundingBox::Data m_split_plane;
};

template<class BoundingBox>
class Partition<BoundingBox, UPPER> {
  public:
    Partition(int axis, typename BoundingBox::Data split_plane) : m_axis(axis), m_split_plane(split_plane) {}
    inline bool operator ()(const BoundingBox &a) const {
      return a.upper(m_axis) < m_split_plane;
    }
  private:
    const int m_axis;
    typename BoundingBox::Data m_split_plane;
};

template<class BoundingBox>
class Partition<BoundingBox, MIDDLE> {
  public:
    Partition(int axis, typename BoundingBox::Data split_plane) : m_axis(axis), m_split_plane(split_plane) {}
    inline bool operator ()(const BoundingBox &a) const {
      return a.middle(m_axis) < m_split_plane;
    }
  private:
    const int m_axis;
    typename BoundingBox::Data m_split_plane;
};

template<class BoundingBox>
class Partition<BoundingBox, LENGTH> {
  public:
    Partition(int axis, typename BoundingBox::Data split_plane) : m_axis(axis), m_split_plane(split_plane) {}
    inline bool operator ()(const BoundingBox &a) const {
      return a.length(m_axis) < m_split_plane;
    }
  private:
    const int m_axis;
    typename BoundingBox::Data m_split_plane;
};

template<class BoundingBox>
class Partition<BoundingBox, KEY> {
  public:
    Partition(typename BoundingBox::Key split_plane) : m_split_plane(split_plane) {}
    inline bool operator ()(const BoundingBox &a) const {
      return a.key < m_split_plane;
    }
  private:
    typename BoundingBox::Key m_split_plane;
};

} // namespace compare
} // namespace box
} // namespace search
} // namespace stk

#endif // stk_search_BoundingBoxCompare_hpp
