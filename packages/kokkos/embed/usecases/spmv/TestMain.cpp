
namespace Test {

int test_host();
int test_cuda();

}

int main()
{
  Test::test_host();

#if defined(HAVE_CUDA)
  Test::test_cuda();
#endif

  return 0 ;
}

