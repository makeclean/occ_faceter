# OCC Faceter

OCC Faceter is a coolection of tools used for the production, modification and serialisation of geometries for [DAGMC](https://svalinn.github.io/DAGMC/)


### Installation

OCC Faceter requires [CGAL](https://cgal.org/), [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) and [OCE](https://github.com/tpaviot/oce) or [OCC](https://www.opencascade.com) to run

Install the dependencies and start using the tool.

```sh
$ sudo apt-get install libcgal-dev
$ sudo apt-get install libocc*-dev occ*
$ sudo apt-get install liblapack-dev libhdf5-dev
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=..
$ make
$ make install
```

### Usage

After installation the bin directory should contain occ_faceter and other tools. To use occ_faceter you will need a CAD geometry saved in Step file format. Th example below assumes you have a Step file called test.step. The faceting tolerance is set as 1.e-3 and the output file is called test.h5m. As an optional step mbconvert is used to convert the new h5m file into an stl file which can be visualized in paraview or visit  

```
$ cd ../bin
$ ./occ_faceter test.step 1.e-3 test.h5m
$ mbconvert test.h5m test.stl
```



