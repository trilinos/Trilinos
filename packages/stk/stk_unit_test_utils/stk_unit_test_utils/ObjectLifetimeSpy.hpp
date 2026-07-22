#ifndef STK_OBJECTLIFETIMESPY_HPP
#define STK_OBJECTLIFETIMESPY_HPP

namespace stk::unit_test_util {

int objectLifetimeSpy_getNumConstructions();
int objectLifetimeSpy_getNumCopyConstructions();
int objectLifetimeSpy_getNumMoveConstructions();
int objectLifetimeSpy_getNumDestructions();
int objectLifetimeSpy_getNumCopyAssignments();
int objectLifetimeSpy_getNumMoveAssignments();

void objectLifetimeSpy_clearCounts();
bool objectLifetimeSpy_checkBalancedConstructionsDestructions();

class ObjectLifetimeSpy
{
public:
  ObjectLifetimeSpy();
  explicit ObjectLifetimeSpy(int i);
  virtual ~ObjectLifetimeSpy();

  ObjectLifetimeSpy(const ObjectLifetimeSpy& rhs);
  ObjectLifetimeSpy(ObjectLifetimeSpy&& rhs) noexcept;
  ObjectLifetimeSpy& operator=(const ObjectLifetimeSpy& rhs);
  ObjectLifetimeSpy& operator=(ObjectLifetimeSpy&& rhs);

  int value() const { return m_i; }
  void set_value(int newValue) { m_i = newValue; }

private:
  int m_i;
  static constexpr bool m_printEvents = false;  // Change to true for better output to understand event ordering
};

}

#endif // STK_OBJECTLIFETIMESPY_HPP
