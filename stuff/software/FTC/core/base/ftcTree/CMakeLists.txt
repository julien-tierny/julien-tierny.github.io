find_package(Boost REQUIRED)

ttk_add_base_library(ftcTree
   SOURCES FTCTree.cpp FTCTree_CT.cpp FTCTree_MT.cpp Segmentation.cpp
   HEADERS AtomicUF.h AtomicVector.h DataTypes.h FTCCommon.h FTCTree.h
   FTCTree_CT.h FTCTree_MT.h FTCTree_Template.h Node.h Propagation.h
   Segmentation.h SharedData.h SuperArc.h
   LINK triangulation geometry)
target_include_directories(ftcTree PUBLIC ${Boost_INCLUDE_DIR})

option(TTK_ENABLE_FTC_TREE_PROCESS_SPEED "Enable FTC tree process speed" OFF)
option(TTK_ENABLE_FTC_TREE_STATS_TIME    "Enable FTC tree stats on tasks processing" OFF)
option(TTK_ENABLE_FTC_EARLY_EXIT         "Enable FTC tree early exit for tests" OFF)
option(TTK_ENABLE_FTC_COUNT_INIT         "Include FTC init time in the total" OFF)
option(TTK_DISABLE_FTC_PRIORITY          "Disable FTC tasks priority" ON)
option(TTK_DISABLE_FTC_EARLY_TRUNK       "Disable FTC early trunk detection" OFF)
option(TTK_DISABLE_FTC_MIX_MT            "Disable FTC mixed join/split trees" OFF)
option(TTK_DISABLE_FTC_PARA_COMBINE      "Disable FTC new parallel combine" OFF)
option(TTK_DISABLE_FTC_TRUNK_COMBINE     "Disable FTC trunk process during combine" OFF)
mark_as_advanced(TTK_ENABLE_FTC_TREE_STATS_TIME TTK_ENABLE_FTC_EARLY_EXIT
   TTK_ENABLE_FTC_COUNT_INIT TTK_DISABLE_FTC_PRIORITY
   TTK_DISABLE_FTC_EARLY_TRUNK TTK_DISABLE_FTC_MIX_MT
   TTK_DISABLE_FTC_PARA_COMBINE TTK_DISABLE_FTC_TRUNK_COMBINE)

if (TTK_ENABLE_FTC_TREE_PROCESS_SPEED)
   target_compile_definitions(ftcTree PUBLIC TTK_ENABLE_FTC_TREE_PROCESS_SPEED)
endif()

if (TTK_ENABLE_FTC_TREE_STATS_TIME)
   target_compile_definitions(ftcTree PUBLIC TTK_ENABLE_FTC_TREE_STATS_TIME)
endif()

if (TTK_ENABLE_FTC_COUNT_INIT)
   target_compile_definitions(ftcTree PUBLIC TTK_ENABLE_FTC_COUNT_INIT)
endif()

if (TTK_ENABLE_FTC_EARLY_EXIT)
   target_compile_definitions(ftcTree PUBLIC TTK_ENABLE_FTC_EARLY_EXIT)
endif()

if (TTK_DISABLE_FTC_PRIORITY)
   target_compile_definitions(ftcTree PUBLIC TTK_DISABLE_FTC_PRIORITY)
endif()

if (TTK_DISABLE_FTC_EARLY_TRUNK)
   target_compile_definitions(ftcTree PUBLIC TTK_DISABLE_FTC_EARLY_TRUNK)
endif()

if (TTK_DISABLE_FTC_MIX_MT)
   target_compile_definitions(ftcTree PUBLIC TTK_DISABLE_FTC_MIX_MT)
endif()

if (TTK_DISABLE_FTC_PARA_COMBINE)
   target_compile_definitions(ftcTree PUBLIC TTK_DISABLE_FTC_PARA_COMBINE)
endif()

if (TTK_DISABLE_FTC_TRUNK_COMBINE)
   target_compile_definitions(ftcTree PUBLIC TTK_DISABLE_FTC_TRUNK_COMBINE)
endif()

# TODO : remove this when tests are done
set(TTK_FTC_GRAINSIZE_SEGM_ARCS 50 CACHE INTEGER "buildSegmentation arc size")
set(TTK_FTC_GRAINSIZE_SEGM_VERTS 400000 CACHE INTEGER "buildSegmentation vertices chunk size")
set(TTK_FTC_GRAINSIZE_TRUNK_ARCS 650 CACHE INTEGER "trunk arcs close")
set(TTK_FTC_GRAINSIZE_TRUNK_VERTS 50000 CACHE INTEGER "trunk vertices process")
set(TTK_FTC_GRAINSIZE_COMBINE_MAIN 10000 CACHE INTEGER "Combine grainSize in vertices")

mark_as_advanced(TTK_FTC_GRAINSIZE_SEGM_ARCS TTK_FTC_GRAINSIZE_SEGM_VERTS TTK_FTC_GRAINSIZE_TRUNK_ARCS
   TTK_FTC_GRAINSIZE_TRUNK_VERTS TTK_FTC_GRAINSIZE_COMBINE_MAIN)

target_compile_definitions(ftcTree PUBLIC
   TTK_FTC_GRAINSIZE_SEGM_ARCS=${TTK_FTC_GRAINSIZE_SEGM_ARCS}
   TTK_FTC_GRAINSIZE_SEGM_VERTS=${TTK_FTC_GRAINSIZE_SEGM_VERTS}
   TTK_FTC_GRAINSIZE_TRUNK_ARCS=${TTK_FTC_GRAINSIZE_TRUNK_ARCS}
   TTK_FTC_GRAINSIZE_TRUNK_VERTS=${TTK_FTC_GRAINSIZE_TRUNK_VERTS}
   TTK_FTC_GRAINSIZE_COMBINE_MAIN=${TTK_FTC_GRAINSIZE_COMBINE_MAIN}
   )

