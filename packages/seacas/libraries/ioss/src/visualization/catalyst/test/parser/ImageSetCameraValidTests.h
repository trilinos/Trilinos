// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __IMAGE_SET_CAMERA_VALID_TESTS_H
#define __IMAGE_SET_CAMERA_VALID_TESTS_H

#include <string>

std::string ImagesetCameraLookAtAbsolutePoint = R"(
    begin catalyst
      begin imageset fooImageset
        look at absolute point = 1.1 2 3e-8
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraLookAtRelativePoint = R"(
     begin catalyst
      begin imageset fooImageset
        look at relative point = 1.5 0.5 0.5
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraLookAtElement = R"(
     begin catalyst
      begin imageset fooImageset
        look at element = 17
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraLookAtNode = R"(
     begin catalyst
      begin imageset fooImageset
        look at node = 20
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraLookAtRelativeDistance = R"(
     begin catalyst
      begin imageset fooImageset
        look at relative distance = 2.0
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraLookAtAbsoluteDistance = R"(
     begin catalyst
      begin imageset fooImageset
        look at absolute distance = 15.0
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraUpVector = R"(
     begin catalyst
      begin imageset fooImageset
        up vector = 0 1 2
        image size = 800 450
      end imageset
    end
    )";

std::string ImagesetCameraFOV = R"(
     begin catalyst
      begin imageset fooImageset
        camera fov = 45
        look at relative point = 0.0 0.0 0.0
        look at absolute distance = 10.0
        image size = 800 450
      end imageset
    end
    )";

#endif /* __IMAGE_SET_CAMERA_VALID_TESTS_H */
