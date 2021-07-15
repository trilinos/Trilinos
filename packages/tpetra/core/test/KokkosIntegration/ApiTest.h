#ifndef __API_TEST_H__
#define __API_TEST_H__

class ApiTest {
 public:
  static ApiTest *getInstance();
  void finalizeInstance();
  int setExpectation(std::string, int);
  int setExpectations(std::map<std::string, int> &exp);
  bool testExpectations();
  void map_zero();
  void pincr(std::string key, int value);
  void incr(std::string key);
  void printAll();

 private:
  ApiTest();
  ~ApiTest();
  std::map<std::string, std::pair<int, int> > counter;
  std::list<std::string> funcs {
      "cudaDeviceSynchronize",
      "cudaMemcpy2DAsync",
      "cudaMemcpy3DAsync",
      "cudaMemcpyAsync",
      "cudaMemcpy",
      "cudaMemcpy2D",
      "cudaMemcpy2DArrayToArray",
      "cudaMemcpy2DFromArray",
      "cudaMemcpy2DFromArrayAsync",
      "cudaMemcpy2DToArray",
      "cudaMemcpy2DToArrayAsync"
      "cudaMemcpy3D",
      "cudaMemcpy3DPeer",
      "cudaMemcpy3DPeerAsync",
      "cudaMemcpyFromSymbol"
      "cudaMemcpyFromSymbolAsync",
      "cudaMemcpyPeer",
      "cudaMemcpyPeerAsync",
      "cudaMemcpyToSymbol",
      "cudaMemcpyToSymbolAsync" };
};

#endif
