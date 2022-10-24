// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "CatalystTestFixture.h"
#include "catch.hpp"

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidCamera", "[imageSetInvalid]")
{

  std::string ps = R"(
     begin catalyst
      begin imageset foo
        camera = doesntexist
      end imageset
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidRepresentation", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        representation = doesntexist
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidOperation", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        operation = doesntexist
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidImageBasename", "[notWorking]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image basename = "$%/name"
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidBaseDirectory", "[notWorking]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image base directory = "$%/name"
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidBaseDirectory2", "[notWorking]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image base directory = "/"
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidBaseDirectory3", "[notWorking]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image base directory = "/foo"
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetTwoImageFormats", "[notWorking]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image format = jpg
        image format = png
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidImageFormat", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image format = foo
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidImageSize1", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image size = 10000000 10000000
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidImageSize2", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset foo
        image size = 1000 0
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidOperationShortcut", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin clip clipNone
        relative point on plane = -0.5 0 0
        plane normal = 1 0 0
        side to keep = positive
      end operation
      begin imageset foo
        slice = point 1 2 3 normal 4 5 6
        operation = clipNone
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidCameraShortcut", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin camera fooCamera
      end camera
      begin imageset foo
        look at absolute point = 1 2 3
        camera = fooCamera
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidRepresentationShortcut", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin representation fooRep
      end representation
      begin imageset foo
        representation = fooRep
        color by scalar = mass
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidCameraLookDirection", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset fooImageset
        look direction = 1 2 3
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidCameraShortcut1", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset fooImageset
        camera at absolute point = 0 0 0
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}

TEST_CASE_METHOD(CatalystTestFixture, "ImagesetInvalidCameraShortcut2", "[imageSetInvalid]")
{

  std::string ps = R"(
    begin catalyst
      begin imageset fooImageset
        image name addon = foo
      end imageset
    end
     )";

  checkPhactoriStringInvalidParse(ps);
}
