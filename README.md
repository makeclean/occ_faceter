# OCC Faceter

OCC Faceter is a coolection of tools used for the production, modification and serialisation of geometries for [DAGMC](https://svalinn.github.io/DAGMC/)


### Installation

OCC Faceter requires [CGAL](https://cgal.org/), [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) and [OCE](https://github.com/tpaviot/oce) or [OCC](https://www.opencascade.com) to run

Install the dependencies and start using the tool.

```sh
$ sudo apt-get install libcgal-dev
$ sudo apt-get install libocc*-dev
$ sudo apt-get install occ*
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=..
$ make
$ make installl
```
