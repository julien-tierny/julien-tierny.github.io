The Topology ToolKit (TTK)
=
https://topology-tool-kit.github.io/


# DATA DISCLAIMER
-----------------

The data sets contained in this archive are the property of their copyright
owners. They are made available in this archive for research purpose only. They
should not be redistributed publicly without their owner's permission.

# INSTRUCTIONS:
---------------

The Fibonacci Task-based Reeb graph (FTR) algorithm has been implemented in
the open-source platform
[TTK](https://topology-tool-kit.github.io/index.html).
Installation instructions can be found
[here](https://topology-tool-kit.github.io/installation.html), however notice
that you will require to get the latest version of TTK directly for the git
repository ([here](https://github.com/topology-tool-kit/ttk)) and not from the
tarball in order to get the FTR filter.

If you are an advanced user and you do not wish to install ParaView (TTK's
main user interface), you still have the possibility to install TTK without
ParaView support (although, you will need to have VTK installed). In this
case, you can skip steps 3, 4 and 5 of the installation procedure.
In order to configure TTK to be build only against VTK and not Paraview, you
need to give `cmake` additional arguments:
```bash
$ cmake-gui .. -DTTK_BUILD_PARAVIEW_PLUGINS=OFF
```
If your VTK is not installed system-wide, you can set `-DVTK_DIR` to point to
the build folder of your VTK library.

At the end of the installation process, you should have a *ttkFTRGraphCmd*
executable file in the *build/standalone/FTRGraph/cmd/* folder.

# REPRODUCING THE RESULTS FROM THE PAPER:
-----------------------------------------

The archive *data.tgz* contains two data sets, which are two less sub-sampled
versions of those used in the paper:
   * *starKnot.vtu* (2D)
   * *elephant.vtu* (3D).

You can extract those datasets in the directory containing the FTR executable: (*build/standalone/FTRGraph/cmd*)
Then:

```
$ ./ttkFTRGraphCmd -i starKnot.vtu
[Common]  _____ _____ _  __                       __  __    ____   ___  _  ___
[Common] |_   _|_   _| |/ /                      / /__\ \  |___ \ / _ \/ |/ _ \
[Common]   | |   | | | ' /                      | |/ __| |   __) | | | | | (_) |
[Common]   | |   | | | . \                      | | (__| |  / __/| |_| | |\__, |
[Common]   |_|   |_| |_|\_\                     | |\___| | |_____|\___/|_|  /_/
[Common]                                         \_\  /_/
[Common] Welcome!
[ttkProgramBase] Reading input data...
[ttkProgramBase]   done! (read 2400000 vertices, 4800000 cells)
[OneSkeleton] Edge-list built in 0.434075 s. (7200000 edges, 1 thread(s)).
[ZeroSkeleton] One-skeleton built in 1.02769 s. (12 thread(s)).
[ZeroSkeleton] Vertex edges built in 0.420243 s. (1 thread(s)).
[ZeroSkeleton] Vertex stars built in 0.230016 s. (1 thread(s)).
[ThreeSkeleton] Cell edges built in 0.258343 s. (12 thread(s)).
[ttkFTRGraph] Starting computation on field: Elevation
[FTR Graph]: thread number: 12
[FTR Graph]: debug lvl: 3
[FTR Graph]: segmentation: true
[FTR Graph]: sampling level: 0
[FTR Graph]: alloc time: 0.338025
[FTR Graph]: init time: 0.257757
[FTR Graph]: sort time: 0.035732
[FTR Graph]: simplices sort time: 0.213756
[FTR Graph]: 1434 leaves
[FTR Graph]: leaf search time 0.259461
[FTR Graph]: sweepFrowSeeds time: 1.67393
[FTR Graph]: build time: 1.93382
[FTR Graph]: *TOTAL* time: 2.18823
[ttkFTRGraph] Memory usage: 6226.56 MB.
[ttkProgramBase] Saving output file `output_port#0.vtu'...
[ttkProgramBase] Saving output file `output_port#1.vtu'...
[ttkProgramBase] Saving output file `output_port#2.vtu'...
[Common] Goodbye :)
```

# Program options and arguments:
--------------------------------

From the directory *build/standalone/FTRGraph/cmd/*, running the following command line will
display the program options and arguments:

```
 $./ttkFTRGraphCmd
[CommandLine] Missing mandatory argument:
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]
[CommandLine] Usage:
[CommandLine]   ./ttkFTRGraphCmd
[CommandLine] Argument(s):
[CommandLine]   [-d <Global debug level (default: 3)>]
[CommandLine]   [-t <Global thread number (default: 12)>]
[CommandLine]   [-f <Scalar field id (default: 0)>]
[CommandLine]   -i <{Input data-sets (*.vti, *vtu, *vtp)}>
[CommandLine]   [-o <Output file name base (no extension) (default: `output')>]
[CommandLine] Option(s):
[CommandLine]   [-s: Single sweep (default: 0)]
```

# Using other data-sets:
------------------------

This plugins use VTK Unstructured Grid (vtu) or VTK Image Data (vti) files.
In the case of an image data, an implicit triangulation is built.
(The use of a Reeb graph algorithm is more relevant on unstructured meshes as
regular grids are guaranteed to be simply connected, so contour tree
algorithms may be preferred.)

The scalar field must be a point data, with only one component.
