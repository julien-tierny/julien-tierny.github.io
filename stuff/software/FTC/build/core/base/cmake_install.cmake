# Install script for directory: /home/julien/Pro/git/web/stuff/software/FTC/core/base

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
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk/TTKBaseConfig.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk/TTKBaseConfig.cmake"
         "/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/CMakeFiles/Export/lib/cmake/ttk/TTKBaseConfig.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk/TTKBaseConfig-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk/TTKBaseConfig.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk" TYPE FILE FILES "/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/CMakeFiles/Export/lib/cmake/ttk/TTKBaseConfig.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ttk" TYPE FILE FILES "/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/CMakeFiles/Export/lib/cmake/ttk/TTKBaseConfig-release.cmake")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/abstractTriangulation/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/common/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/explicitTriangulation/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/ftcTree/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/geometry/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/implicitTriangulation/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/skeleton/cmake_install.cmake")
  include("/home/julien/Pro/git/web/stuff/software/FTC/build/core/base/triangulation/cmake_install.cmake")

endif()

