The Topology ToolKit (TTK)
=
https://topology-tool-kit.github.io/


# DATA DISCLAIMER
-----------------

The data sets contained in this archive are the property of their copyright
owners. They are made available in this archive for research purpose only. They
should not be redistributed publicly without their owner's permission.

# CODE DISCLAIMER
-----------------

This code is provided "as is" without warranty of any kind. It is made available
under the licensing terms described in the file 'LICENSE'.

# SOFTWARE REQUIREMENTS:
------------------------

This code requires:

* the VTK development libraries to be installed on the system
  (version 6.1 or higher). Please see http://www.vtk.org/ for installation
  instructions.

* A recent version of CMake (http://www.cmake.org) should be installed (version
  2.8 or higher).

* A recent version of Boost (http://www.boost.org/) should be installed (version
  1.49.0 or higher).

# INSTRUCTIONS:
---------------

The following instructions are for Unix environments (MacOS, Linux, etc.). For
Windows users, please follow the usual procedure to import a CMake project into
VisualStudio (please be careful though, important directory paths may need to
be changed for the program to run properly).

## Standalone only:

Compiling the code:

After decompressing the archive, enter its root directory.
From there, enter the following commands (omit the "$" character):

``` bash
$ mkdir build
$ cd build
$ cmake --clean-first .. -DCMAKE_INSTALL_PREFIX=.
$ make install
```

If you have built your own version of VTK, you can tell cmake where
to find it by adding: `-DVTK_DIR=<vtk-location>/build`

The standalone is:  *build/standalone/FTCTree/cmd/ftcTreeCmd*

# REPRODUCING THE RESULTS FROM THE PAPER:
-----------------------------------------

From the directory containing the executable (*build/standalone/FTCTree/cmd*)

```
$ ./ftcTreeCmd -g ../../../../data.vti -f 1 -d 3

[Common]  _____ _____ _  __                       __  __    ____   ___  _  ___
[Common] |_   _|_   _| |/ /                      / /__\ \  |___ \ / _ \/ |( _ )
[Common]   | |   | | | ' /                      | |/ __| |   __) | | | | |/ _ \
[Common]   | |   | | | . \                      | | (__| |  / __/| |_| | | (_) |
[Common]   |_|   |_| |_|\_\                     | |\___| | |_____|\___/|_|\___/
[Common]                                         \_\  /_/
[Common] Welcome!
[Editor] Reading input mesh...
[Editor]   done! (read 2097152 vertices, 2048383 cells)
[ImplicitTriangulation] The getVertex*() requests are accelerated.
Launch on field : ftle
------------
[FTC] number of threads : 12
* debug lvl  : 3
* tree type  : Contour
------------
[FTC] alloc in           0.0833402
[FTC] init in            0.0330589
[FTC] sort step in       0.0469201
[FTC] leafSearch in      0.0183449
[FTC] leaf Growth JT in  0.0302839
[FTC] leaf Growth ST in  0.035357
[FTC] trunk JT in        0.0201309
[FTC] trunk ST in        0.0310669
[FTC] segment JT in      0.0211661
[FTC] segment ST in      0.016603
[FTC] merge trees  in    0.083317
[FTC] parallel combine in 0.0144999
[FTC] trunk combine in   0.0142581
[FTC] Combine end in     0.028811
[FTC] combine full in    0.050549
[FTC] build tree in      0.152242
[FTC] Total  in          0.199186
[ttkFTCTree] Memory usage: 7852.91 MB.
[Common] Goodbye :)
```

# Program options and arguments:
--------------------------------

From the directory *build/standalone/FTCTree/cmd/*, running the following command line will
display the program options and arguments:

```
./ftcTreeCmd
[Common]  _____ _____ _  __                       __  __    ____   ___  _  ___
[Common] |_   _|_   _| |/ /                      / /__\ \  |___ \ / _ \/ |( _ )
[Common]   | |   | | | ' /                      | |/ __| |   __) | | | | |/ _ \
[Common]   | |   | | | . \                      | | (__| |  / __/| |_| | | (_) |
[Common]   |_|   |_| |_|\_\                     | |\___| | |_____|\___/|_|\___/
[Common]                                         \_\  /_/
[Common] Welcome!
[CommandLine] Missing mandatory argument:
[CommandLine]   -g <Path to the input 3D grid (default: `')>
[CommandLine]
[CommandLine] Usage:
[CommandLine]   ./ftcTreeCmd
[CommandLine] Argument(s):
[CommandLine]   [-d <Global debug level (default: 3)>]
[CommandLine]   [-t <Thread number (default: 12)>]
[CommandLine]   -g <Path to the input 3D grid (default: `')>
[CommandLine]   [-f <Field identifier (default: 0)>]
[CommandLine]   [-T <type of tree : 2 is CT (default: 2)>]
[CommandLine] Option(s):
```

# Using other data-sets:
------------------------

This plugins use VTK Unstructured Grid (vtu) or VTK Image Data (vti) files.
In the case of an image data, an implicit triangulation is built.

The scalar field must be a point data, with only one component of type float.


# Using this implementation from your own code:
-----------------------------------------------

Please see the source files *standalone/FTCTree/cmd/Editor.cpp* for usage examples.
