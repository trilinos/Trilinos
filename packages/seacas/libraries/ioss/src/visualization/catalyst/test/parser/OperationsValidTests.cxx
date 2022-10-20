// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "OperationsValidTests.h"
#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture, "OperationBaseline", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  checkPhactoriStringValidParse(OperationBaseline, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClip", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]                   = "clip";
  fooOperation["input"]                  = "clipNone";
  dj["operation blocks"]["fooOperation"] = fooOperation;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClip, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipAbsolutePoint", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]  = "clip";
  fooOperation["input"] = "clipNone";
  Json::Value plane;
  plane[0]                                = 3.0;
  plane[1]                                = 0.1;
  plane[2]                                = 0.2;
  fooOperation["absolute point on plane"] = plane;
  Json::Value normal;
  normal[0]                    = 1.0;
  normal[1]                    = 1.0;
  normal[2]                    = 0.0;
  fooOperation["plane normal"] = normal;

  dj["operation blocks"]["fooOperation"] = fooOperation;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClipAbsolutePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipRelativePoint", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]  = "clip";
  fooOperation["input"] = "clipNone";
  Json::Value plane;
  plane[0]                                = 0.25;
  plane[1]                                = 0.3;
  plane[2]                                = 0.5;
  fooOperation["relative point on plane"] = plane;
  Json::Value normal;
  normal[0]                              = 1.0;
  normal[1]                              = 1.0;
  normal[2]                              = 0.0;
  fooOperation["plane normal"]           = normal;
  dj["operation blocks"]["fooOperation"] = fooOperation;

  Json::Value &camera                 = dj["camera blocks"]["singletestcamera"];
  camera["look at absolute distance"] = 10.0;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClipRelativePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipNode", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]          = "clip";
  fooOperation["input"]         = "clipNone";
  fooOperation["node on plane"] = 20;
  Json::Value normal;
  normal[0]                              = 1.0;
  normal[1]                              = 1.0;
  normal[2]                              = 0.0;
  fooOperation["plane normal"]           = normal;
  dj["operation blocks"]["fooOperation"] = fooOperation;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClipNode, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipElement", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]             = "clip";
  fooOperation["input"]            = "clipNone";
  fooOperation["element on plane"] = 17;
  Json::Value normal;
  normal[0]                              = 1.0;
  normal[1]                              = 0.25;
  normal[2]                              = 0.0;
  fooOperation["plane normal"]           = normal;
  dj["operation blocks"]["fooOperation"] = fooOperation;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClipElement, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipPlaneNormal", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  Json::Value fooOperation;
  fooOperation["type"]  = "clip";
  fooOperation["input"] = "clipNone";
  Json::Value normal;
  normal[0]                    = 1.0;
  normal[1]                    = 2.0;
  normal[2]                    = 3.0;
  fooOperation["plane normal"] = normal;
  Json::Value point;
  point[0]                                = 3.0;
  point[1]                                = 0.1;
  point[2]                                = 0.2;
  fooOperation["absolute point on plane"] = point;
  dj["operation blocks"]["fooOperation"]  = fooOperation;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  imageset["operation"] = "fooOperation";
  checkPhactoriStringValidParse(OperationClipPlaneNormal, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "OperationClipPlaneThreePoints1", "[operationsValid]")
{
  Json::Value dj = getDefaultOperationsJSON();
  dj["operation blocks"].removeMember("clipNone");
  Json::Value eb;
  eb["type"] = "extractblock";
  Json::Value ib;
  ib[0]                         = "block_1";
  eb["include blocks"]          = ib;
  dj["operation blocks"]["eb1"] = eb;

  Json::Value fooOperation;
  fooOperation["type"]                = "clip";
  fooOperation["input"]               = "eb1";
  fooOperation["plane specification"] = "three points";
  fooOperation["side to keep"]        = "negative";
  Json::Value planeA;
  planeA[0]                             = "max";
  planeA[1]                             = "scalar";
  planeA[2]                             = "VON_MISES";
  planeA[3]                             = "eb1";
  fooOperation["data point on plane A"] = planeA;
  Json::Value planeB;
  planeB[0]                             = "min";
  planeB[1]                             = "scalar";
  planeB[2]                             = "VON_MISES";
  planeB[3]                             = "eb1";
  fooOperation["data point on plane B"] = planeB;
  Json::Value planeC;
  planeC[0]                                 = 0.0;
  planeC[1]                                 = 10.0;
  planeC[2]                                 = 0.0;
  fooOperation["relative point on plane C"] = planeC;
  dj["operation blocks"]["fooOperation"]    = fooOperation;

  Json::Value clipop2;
  clipop2["type"]                = "clip";
  clipop2["plane specification"] = "three points";
  Json::Value planeA2;
  planeA2[0]                       = "max";
  planeA2[1]                       = "scalar";
  planeA2[2]                       = "VON_MISES";
  planeA2[3]                       = "eb1";
  clipop2["data point on plane A"] = planeA2;
  Json::Value planeB2;
  planeB2[0]                           = 0.0;
  planeB2[1]                           = 0.0;
  planeB2[2]                           = 10.0;
  clipop2["relative point on plane B"] = planeB2;
  Json::Value planeC2;
  planeC2[0]                           = 0.0;
  planeC2[1]                           = 10.0;
  planeC2[2]                           = 0.0;
  clipop2["relative point on plane C"] = planeC2;
  dj["operation blocks"]["clipop2"]    = clipop2;

  Json::Value &camera = dj["camera blocks"]["singletestcamera"];
  camera.removeMember("look direction");
  Json::Value point1;
  point1[0]                          = 15.0;
  point1[1]                          = 10.0;
  point1[2]                          = 10.0;
  camera["camera at absolute point"] = point1;
  Json::Value point2;
  point2[0]                        = 5.0;
  point2[1]                        = 0.0;
  point2[2]                        = 0.0;
  camera["look at absolute point"] = point2;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  is2(imageset);
  Json::Value  is3(imageset);
  imageset["operation"]       = "fooOperation";
  imageset["color by scalar"] = "VON_MISES";
  imageset["show edges"]      = true;
  is2["operation"]            = "clipop2";
  is2["show edges"]           = true;
  is2["color by scalar"]      = "VON_MISES";
  is3["show edges"]           = true;
  is3["color by scalar"]      = "VON_MISES";
  is3.removeMember("operation");

  dj["imageset blocks"]["is2"] = is2;
  dj["imageset blocks"]["is3"] = is3;

  checkPhactoriStringValidParse(OperationClipPlaneThreePoints1, dj);
}
