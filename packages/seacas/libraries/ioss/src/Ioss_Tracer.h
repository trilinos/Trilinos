#ifndef IOSS_Tracer_h
#define IOSS_Tracer_h
namespace Ioss {
  class Tracer
  {
  public:
    Tracer(const char *function);
    ~Tracer();

  private:
    static int level;
  };
}
#endif
