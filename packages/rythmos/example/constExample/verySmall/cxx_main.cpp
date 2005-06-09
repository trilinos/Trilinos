
class Bar
{
  public:
    Bar() {};
    ~Bar() {};
  protected:
    double t_;
};

class Foo
{
  public:
    Foo() { };
    ~Foo() { };
    const Bar &getBar() { return(B_); };
  protected:
    Bar B_;
};

const Bar &somefunction(const Foo &F)
{
  const Bar &B = F.getBar(); // This doesn't work.  Why?
  return(B);
};

int main(int argc, char *argv[])
{
  Foo F;
  const Bar &B = somefunction(F);

  return(0);
}
