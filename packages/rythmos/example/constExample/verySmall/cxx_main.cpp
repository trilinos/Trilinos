
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
//    const Bar &getBar() const { return(B_); }; // this allows the line below to work.
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
  const Bar &B1 = F.getBar(); // This works.
  const Bar &B2 = somefunction(F);

  return(0);
}
