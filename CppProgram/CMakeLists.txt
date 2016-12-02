# Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
# Authors: Ugo Pattacini
# CopyPolicy: Released under the terms of the GNU GPL v2.0.

cmake_minimum_required(VERSION 2.8.9)
project(learningTrajectories)

find_package(YARP         REQUIRED)
find_package(ICUBcontrib  REQUIRED)
find_package(hapticdevice REQUIRED)
find_package(gazebo REQUIRED)

list(APPEND CMAKE_MODULE_PATH ${YARP_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${ICUBCONTRIB_MODULE_PATH})

include(ICUBcontribOptions)
include(ICUBcontribHelpers)
icubcontrib_set_default_prefix()


include_directories(${YARP_INCLUDE_DIRS} ${hapticdevice_INCLUDE_DIRS})

add_library(fake_geomagic_driver fakeGeomagicDriver.h fakeGeomagicDriver.cpp)
target_link_libraries(fake_geomagic_driver ${YARP_LIBRARIES}  )
    
add_definitions(-D_USE_MATH_DEFINES)


include_directories(${GAZEBO_INCLUDE_DIRS})
link_directories(${GAZEBO_LIBRARY_DIRS})
list(APPEND CMAKE_CXX_FLAGS "${GAZEBO_CXX_FLAGS}")

#add_executable(test_simple test_simple.cpp cartesianClient.cpp)
#target_link_libraries(test_simple ${YARP_LIBRARIES} pthread fake_geomagic_driver)


add_executable(recordTrajectoriesWithGeomagic recordTrajectoriesWithGeomagic.cpp cartesianClient.cpp)
target_link_libraries(recordTrajectoriesWithGeomagic ${YARP_LIBRARIES} pthread fake_geomagic_driver)


add_executable(replayTrajectories replayTrajectories.cpp cartesianClient.cpp)
target_link_libraries(replayTrajectories ${YARP_LIBRARIES} pthread)

add_executable(recordTrajectoriesWithiCub recordTrajectoriesWithiCub.cpp cartesianClient.cpp)
target_link_libraries(recordTrajectoriesWithiCub ${YARP_LIBRARIES} pthread SDL )

#add_executable(recordTrajectoriesWithiCub recordTrajectoriesWithiCubV2.cpp)
#target_link_libraries(recordTrajectoriesWithiCub ${YARP_LIBRARIES} pthread SDL )


install(TARGETS recordTrajectoriesWithGeomagic DESTINATION bin)
install(TARGETS replayTrajectories DESTINATION bin)
#install(TARGETS test_simple DESTINATION bin)
install(TARGETS recordTrajectoriesWithiCub DESTINATION bin)
