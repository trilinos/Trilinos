// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"
#include "ImageSetOperationValidTests.h"

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold1", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "scalar";
  th[1] = "VON_MISES";
  th[2] = "keep between";
  Json::Value kb;
  kb[0]                 = 28.5;
  kb[1]                 = 30.5;
  th[3]                 = kb;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold2", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "vector magnitude";
  th[1] = "displ";
  th[2] = "keep between";
  Json::Value kb;
  kb[0]                 = 0.003;
  kb[1]                 = 0.008;
  th[3]                 = kb;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold2, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold3", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "vector component";
  th[1] = "displ_x";
  th[2] = "keep between";
  Json::Value kb;
  kb[0]                 = 0.003;
  kb[1]                 = 0.008;
  th[3]                 = kb;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold3, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold4", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "tensor component";
  th[1] = "stress_xy";
  th[2] = "keep between";
  Json::Value kb;
  kb[0]                 = 1.0;
  kb[1]                 = 1.7;
  th[3]                 = kb;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold4, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold5", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "scalar";
  th[1] = "VON_MISES";
  th[2] = "keep above";
  Json::Value ka;
  ka[0]                 = 28.5;
  th[3]                 = ka;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold5, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationThreshold6", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  th;
  th[0] = "scalar";
  th[1] = "VON_MISES";
  th[2] = "keep below";
  Json::Value kb;
  kb[0]                 = 28.5;
  th[3]                 = kb;
  imageset["threshold"] = th;

  checkPhactoriStringValidParse(ImagesetOperationThreshold6, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationBoxClip1", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  bc;
  Json::Value  center;
  center[0] = 4.0;
  center[1] = 1.0;
  center[2] = -0.61;
  bc[0]     = center;
  Json::Value extents;
  extents[0]          = 6.0;
  extents[1]          = 1.5;
  extents[2]          = 2.0;
  bc[1]               = extents;
  bc[2]               = "keep inside";
  imageset["boxclip"] = bc;

  checkPhactoriStringValidParse(ImagesetOperationBoxClip1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationClip1", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  Json::Value  point;
  point[0] = 0.0;
  point[1] = 1.0;
  point[2] = 2.0;
  c[0]     = point;
  Json::Value normal;
  normal[0]        = 3.0;
  normal[1]        = 4.0;
  normal[2]        = 5.0;
  c[1]             = normal;
  imageset["clip"] = c;

  checkPhactoriStringValidParse(ImagesetOperationClip1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationSlice1", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  s;
  Json::Value  point;
  point[0] = 0.0;
  point[1] = 1.0;
  point[2] = 2.0;
  s[0]     = point;
  Json::Value normal;
  normal[0]         = 3.0;
  normal[1]         = 4.0;
  normal[2]         = 5.0;
  s[1]              = normal;
  imageset["slice"] = s;

  checkPhactoriStringValidParse(ImagesetOperationSlice1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour1", "[ImageSetOperationValid]")
{
  Json::Value  dj     = getDefaultImageSetWithCameraJSON();
  Json::Value &camera = dj["camera blocks"]["singletestcamera"];
  Json::Value  direction;
  direction[0]             = -1.0;
  direction[1]             = -1.0;
  direction[2]             = -1.0;
  camera["look direction"] = direction;
  Json::Value point;
  point[0]                            = 5.0;
  point[1]                            = 1.0;
  point[2]                            = 0.0;
  camera["look at absolute point"]    = point;
  camera["look at absolute distance"] = 10.0;

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "scalar";
  c[1] = "VON_MISES";
  Json::Value valueList;
  valueList[0]        = 30.0;
  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour1, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour2", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "vector magnitude";
  c[1] = "displ";
  Json::Value valueList;
  valueList[0]        = 0.0001;
  valueList[1]        = 0.0021;
  valueList[2]        = 0.0041;
  valueList[3]        = 0.0061;
  valueList[4]        = 0.0081;
  valueList[5]        = 0.011;
  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour2, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour3", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "vector component";
  c[1] = "displ_x";
  Json::Value valueList;
  valueList[0]        = 0.0001;
  valueList[1]        = 0.0021;
  valueList[2]        = 0.0041;
  valueList[3]        = 0.0061;
  valueList[4]        = 0.0081;
  valueList[5]        = 0.011;
  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour3, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour4", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "tensor component";
  c[1] = "stress_xy";
  Json::Value valueList;
  valueList[0]        = 0.1;
  valueList[1]        = 0.5;
  valueList[2]        = 1.0;
  valueList[3]        = 1.5;
  valueList[4]        = 1.8;
  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour4, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour5", "[ImageSetOperationValid]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "vector component";
  c[1] = "displ_x";
  Json::Value valueList;
  valueList[0]  = 0.0001;
  valueList[1]  = 0.0011;
  valueList[2]  = 0.0021;
  valueList[3]  = 0.0031;
  valueList[4]  = 0.0041;
  valueList[5]  = 0.0051;
  valueList[6]  = 0.0061;
  valueList[7]  = 0.0071;
  valueList[8]  = 0.0081;
  valueList[9]  = 0.0091;
  valueList[10] = 0.0101;

  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour5, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetOperationContour6", "[notWorking]")
{
  Json::Value dj = getDefaultImageSetWithCameraJSON();

  Json::Value &imageset = dj["imageset blocks"]["fooImageset"];
  Json::Value  c;
  c[0] = "vector component";
  c[1] = "displ_x";
  Json::Value valueList;
  valueList[0] = 0.0001;
  valueList[1] = 0.001;
  valueList[2] = 0.03;

  c[2]                = "value list";
  c[3]                = valueList;
  imageset["contour"] = c;

  checkPhactoriStringValidParse(ImagesetOperationContour6, dj);
}
