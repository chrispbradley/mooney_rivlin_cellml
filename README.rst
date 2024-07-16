====================
Mooney Rivlin CellML
====================

Finite elasticity example using the Mooney Rivlin constitutive relationship described in CellML.

Building the example
====================

Instructions on how to configure and build with CMake::

  git clone https://github.com/OpenCMISS-Examples/mooney_rivlin_cellml.git
  cd mooney_rivlin_cellml
  mkdir build
  cd build
  cmake -DOpenCMISS_INSTALL_ROOT=/path/to/opencmiss/install ../.
  make  # cmake --build . will also work here and is much more platform agnostic.

Running the example
===================

To run the example execute the following commands from the build directory::

  ./src/fortran/mooney_rivlin_cellml

Prerequisites
=============

This example is self-contained and requires no further pre-requisties other than Iron.

License
=======

Apache 2.0 License
