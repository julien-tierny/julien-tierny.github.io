# Install script for directory: /home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/julien/Pro/git/web/stuff/software/FTC/build")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/ttk" TYPE STATIC_LIBRARY FILES "/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/ftcTree/libftcTree.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ttk/base" TYPE FILE FILES
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/AtomicUF.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/AtomicVector.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/DataTypes.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/FTCCommon.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/FTCTree.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/FTCTree_CT.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/FTCTree_MT.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/FTCTree_Template.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/Node.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/Propagation.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/Segmentation.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/SharedData.h"
    "/home/julien/Pro/git/web/stuff/software/FTC/core/base/ftcTree/SuperArc.h"
    )
endif()

