#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ttk::vtk::ttkFTCTree" for configuration "Release"
set_property(TARGET ttk::vtk::ttkFTCTree APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttk::vtk::ttkFTCTree PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/ttk/libttkFTCTree.so"
  IMPORTED_SONAME_RELEASE "libttkFTCTree.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS ttk::vtk::ttkFTCTree )
list(APPEND _IMPORT_CHECK_FILES_FOR_ttk::vtk::ttkFTCTree "${_IMPORT_PREFIX}/lib/ttk/libttkFTCTree.so" )

# Import target "ttk::vtk::ttkProgramBase" for configuration "Release"
set_property(TARGET ttk::vtk::ttkProgramBase APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttk::vtk::ttkProgramBase PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/ttk/libttkProgramBase.so"
  IMPORTED_SONAME_RELEASE "libttkProgramBase.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS ttk::vtk::ttkProgramBase )
list(APPEND _IMPORT_CHECK_FILES_FOR_ttk::vtk::ttkProgramBase "${_IMPORT_PREFIX}/lib/ttk/libttkProgramBase.so" )

# Import target "ttk::vtk::ttkTextureMapFromField" for configuration "Release"
set_property(TARGET ttk::vtk::ttkTextureMapFromField APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttk::vtk::ttkTextureMapFromField PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/ttk/libttkTextureMapFromField.so"
  IMPORTED_SONAME_RELEASE "libttkTextureMapFromField.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS ttk::vtk::ttkTextureMapFromField )
list(APPEND _IMPORT_CHECK_FILES_FOR_ttk::vtk::ttkTextureMapFromField "${_IMPORT_PREFIX}/lib/ttk/libttkTextureMapFromField.so" )

# Import target "ttk::vtk::ttkTriangulation" for configuration "Release"
set_property(TARGET ttk::vtk::ttkTriangulation APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttk::vtk::ttkTriangulation PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/ttk/libttkTriangulation.so"
  IMPORTED_SONAME_RELEASE "libttkTriangulation.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS ttk::vtk::ttkTriangulation )
list(APPEND _IMPORT_CHECK_FILES_FOR_ttk::vtk::ttkTriangulation "${_IMPORT_PREFIX}/lib/ttk/libttkTriangulation.so" )

# Import target "ttk::vtk::ttkUserInterfaceBase" for configuration "Release"
set_property(TARGET ttk::vtk::ttkUserInterfaceBase APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ttk::vtk::ttkUserInterfaceBase PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/ttk/libttkUserInterfaceBase.so"
  IMPORTED_SONAME_RELEASE "libttkUserInterfaceBase.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS ttk::vtk::ttkUserInterfaceBase )
list(APPEND _IMPORT_CHECK_FILES_FOR_ttk::vtk::ttkUserInterfaceBase "${_IMPORT_PREFIX}/lib/ttk/libttkUserInterfaceBase.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
