// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __CAMERA_VALID_TESTS_H
#define __CAMERA_VALID_TESTS_H

#include <string>

std::string Camera = R"(
     begin catalyst
       begin camera fooCamera
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraLookAtAbsolutePoint = R"(
     begin catalyst
       begin camera fooCamera
         look at absolute point = 1.1 2 3e-8
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraLookDirection = R"(
     begin catalyst
       begin camera fooCamera
         look direction = 1 2 3
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraLookAtRelativeDistance = R"(
     begin catalyst
       begin camera fooCamera
         look at relative distance = 2.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraLookAtAbsoluteDistance = R"(
     begin catalyst
       begin camera fooCamera
         look at absolute distance = 15.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtAbsolutePoint = R"(
     begin catalyst
       begin camera fooCamera
         camera at absolute point = -2.0 3.0 30.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtRelativePoint = R"(
     begin catalyst
       begin camera fooCamera
         camera at relative point = -0.5 1.5 20.0
         look direction = 0.1 -0.1 -1.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtElement = R"(
     begin catalyst
       begin camera fooCamera
         camera at element = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtNode = R"(
     begin catalyst
       begin camera fooCamera
         camera at node = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtElementDisplaced = R"(
     begin catalyst
       begin camera fooCamera
         camera at element displaced = 1 -3.0 10.0 20.0
         look at element = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraAtNodeDisplaced = R"(
     begin catalyst
       begin camera fooCamera
         camera at node displaced = 1 -3.0 10.0 20.0
         look at node = 1
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraUpVector = R"(
     begin catalyst
       begin camera fooCamera
         up vector = 0 1 2
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraFOV = R"(
     begin catalyst
       begin camera fooCamera
         camera fov = 45
         look at relative point = 0.0 0.0 0.0
         look at absolute distance = 10.0
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraImageNameAddon = R"(
     begin catalyst
       begin camera fooCamera
         image name addon = "_foo_"
       end camera
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string CameraParallelProjection1 = R"(
     begin catalyst
       begin camera parallel_projection_cam1
         projection type = parallel
         look direction = -5 -1 -1
       end camera
       begin camera perspective_projection_cam1
         projection type = perspective
         look direction = -5 -1 -1
       end camera
       begin imageset parallel_projection_is1
         camera = parallel_projection_cam1
         image size = 800 450
         image basename = parallel_projection_is1.
       end imageset
       begin imageset perspective_projection_is1
         camera = perspective_projection_cam1
         image size = 800 450
         image basename = perspective_projection_is1.
       end imageset
     end
     )";

std::string MultiCamera8 = R"(
     begin catalyst
       begin multicamera8 fooCamera
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtAbsolutePoint = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at absolute point = 1.1 2 3e-8
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtRelativePoint = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at relative point = 0.5 -0.5 0.5
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtElement = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at element = 17
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtNode = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at node = 20
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtRelativeDistance = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at relative distance = 2.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8LookAtAbsoluteDistance = R"(
     begin catalyst
       begin multicamera8 fooCamera
         look at absolute distance = 15.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8UpVector = R"(
     begin catalyst
       begin multicamera8 fooCamera
         up vector = 0 1 2
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8FOV = R"(
     begin catalyst
       begin multicamera8 fooCamera
         camera fov = 45
         look at relative point = 0.0 0.0 0.0
         look at absolute distance = 10.0
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8ImageNameAddon = R"(
     begin catalyst
       begin multicamera8 fooCamera
         image name addon = "_foo_"
       end multicamera8
       begin imageset fooImageset
         camera = fooCamera
         image size = 800 450
       end imageset
     end
     )";

std::string MultiCamera8ParallelProjection1 = R"(
     begin catalyst
       begin multicamera8 parallel_projection_cam1
         projection type = parallel
       end multicamera8
       begin multicamera8 perspective_projection_cam1
         projection type = perspective
       end multicamera8
       begin imageset parallel_projection_is1
         camera = parallel_projection_cam1
         image size = 800 450
         image basename = parallel_projection_is1.
       end imageset
       begin imageset perspective_projection_is1
         camera = perspective_projection_cam1
         image size = 800 450
         image basename = perspective_projection_is1.
       end imageset
     end
     )";

#endif /* __CAMERA_VALID_TESTS_H */
