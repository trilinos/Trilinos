// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"
#include "ImageSetCameraValidTests.h"

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtAbsolutePoint", "[ImageSetCameraValid]")
{
  Json::Value  dj                       = getDefaultImageSetJSON();
  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["look at absolute point"][0] = 1.1;
  imageset["look at absolute point"][1] = 2.0;
  imageset["look at absolute point"][2] = 3e-8;

  checkPhactoriStringValidParse(ImagesetCameraLookAtAbsolutePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtRelativePoint", "[ImageSetCameraValid]")
{
  Json::Value  dj                       = getDefaultImageSetJSON();
  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["look at relative point"][0] = 1.5;
  imageset["look at relative point"][1] = 0.5;
  imageset["look at relative point"][2] = 0.5;

  checkPhactoriStringValidParse(ImagesetCameraLookAtRelativePoint, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtElement", "[ImageSetCameraValid]")
{
  Json::Value  dj             = getDefaultImageSetJSON();
  Json::Value &imageset       = dj["imageset blocks"]["fooImageset"];
  imageset["look at element"] = 17;

  checkPhactoriStringValidParse(ImagesetCameraLookAtElement, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtNode", "[ImageSetCameraValid]")
{
  Json::Value  dj          = getDefaultImageSetJSON();
  Json::Value &imageset    = dj["imageset blocks"]["fooImageset"];
  imageset["look at node"] = 20;

  checkPhactoriStringValidParse(ImagesetCameraLookAtNode, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtRelativeDistance",
                 "[ImageSetCameraValid]")
{
  Json::Value  dj                       = getDefaultImageSetJSON();
  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["look at relative distance"] = 2.0;

  checkPhactoriStringValidParse(ImagesetCameraLookAtRelativeDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraLookAtAbsoluteDistance",
                 "[ImageSetCameraValid]")
{
  Json::Value  dj                       = getDefaultImageSetJSON();
  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["look at absolute distance"] = 15.0;

  checkPhactoriStringValidParse(ImagesetCameraLookAtAbsoluteDistance, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraUpVector", "[ImageSetCameraValid]")
{
  Json::Value  dj          = getDefaultImageSetJSON();
  Json::Value &imageset    = dj["imageset blocks"]["fooImageset"];
  imageset["up vector"][0] = 0.0;
  imageset["up vector"][1] = 1.0;
  imageset["up vector"][2] = 2.0;

  checkPhactoriStringValidParse(ImagesetCameraUpVector, dj);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetCameraFOV", "[ImageSetCameraValid]")
{
  Json::Value  dj                       = getDefaultImageSetJSON();
  Json::Value &imageset                 = dj["imageset blocks"]["fooImageset"];
  imageset["camera fov"]                = 45.0;
  imageset["look at relative point"][0] = 0.0;
  imageset["look at relative point"][1] = 0.0;
  imageset["look at relative point"][2] = 0.0;
  imageset["look at absolute distance"] = 10.0;

  checkPhactoriStringValidParse(ImagesetCameraFOV, dj);
}
