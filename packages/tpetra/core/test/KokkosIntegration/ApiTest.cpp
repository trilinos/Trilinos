
#include <string>
#include <map>
#include <list>
#include "ApiTest.h"

ApiTest *instance = NULL;

ApiTest *ApiTest::getInstance() {
  if (instance == NULL)
    instance = new ApiTest;
  return instance;
}

void ApiTest::finalizeInstance() {
  if (instance != NULL)
    delete instance;
}

ApiTest::ApiTest() { }

ApiTest::~ApiTest() { }

int ApiTest::setExpectation(std::string key, int value) {
  std::map<std::string, std::pair<int, int> >::iterator found = counter.find(key);
  if (found != counter.end()) {
    found->second.second = value;
    return 0;
  }
  return -1;
}

int ApiTest::setExpectations(std::map<std::string, int> &exp) {
  for (std::map<std::string, int>::iterator it = exp.begin();
       it != exp.end(); it++) {
    std::map<std::string, std::pair<int, int> >::iterator found = counter.find(it->first);
    if (found != counter.end()) {
      found->second.second = it->second;
    } else {
      return -1;
    }
  }
  return 0;
}

bool ApiTest::testExpectations() {
  for (std::map<std::string, std::pair<int, int> >::iterator it = counter.begin();
       it != counter.end(); it++) {
    if (it->second.first != it->second.second)
      return false;
  }
  return true;
}

void ApiTest::map_zero() {
  for (std::list<std::string>::iterator it = ApiTest::funcs.begin(); it != ApiTest::funcs.end(); it++) {
    std::map<std::string, std::pair<int, int> >::iterator found = counter.find(*it);
    if (found != counter.end()) {
      found->second = std::make_pair(0, 0);
    } else {
      counter.insert(std::make_pair(*it, std::make_pair(0, 0)));
    }
  }
}

void ApiTest::pincr(std::string key, int val) {
  std::map<std::string, std::pair<int, int> >::iterator it = counter.find(key);
  if (it == counter.end()) {
    counter.insert(std::make_pair(key, std::make_pair(1, 0)));
  } else {
    it->second.first += val;
  }
}

void ApiTest::incr(std::string key) {
  pincr(key, 1);
}

void ApiTest::printAll() {
  fprintf (stderr, "%-30s %10s\t %10s\n","kernel","actual","expected");
  for (std::map<std::string, std::pair<int, int> >::iterator it = counter.begin(); 
       it != counter.end(); it++) {
    if (it->second.first != 0 || it->second.second != 0) 
      fprintf(stderr, "%-30s %10d\t%10d\n", it->first.c_str(), it->second.first, it->second.second);
  }
}
