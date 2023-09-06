# OCC Faceter

OCC Faceter is a coolection of tools used for the production, modification and serialisation of geometries for [DAGMC](https://svalinn.github.io/DAGMC/)


### Installation

*See also the installation in [./Dockerfile](./Dockerfile).*

OCC Faceter has a number of dependencies, including [CGAL](https://cgal.org/), [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) and [OCC](https://www.opencascade.com).

To install MOAB and its dependencies:

```
sudo apt install libhdf5-dev build-essential gfortran autoconf libtool liblapack-dev
git clone https://bitbucket.org/fathomteam/moab.git
cd moab
autoreconf -fi
mkdir bld; cd bld
../configure --enable-optimize --enable-shared --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial --prefix=/opt/moab
make -j4
make check
sudo mkdir /opt/moab
sudo chown $USER /opt/moab
make install
cd ../..
```

To install OCC:

```sh
sudo apt install libocct*-dev occt*
```

To install occ_faceter and remaining dependencies:

```
sudo apt install cmake libcgal-dev libtbb-dev
git clone https://github.com/makeclean/occ_faceter.git
cd occ_faceter
mkdir build
cd build
export LD_LIBRARY_PATH=/opt/moab/lib
# build occ_faceter (to also build other tools, add -DBUILD_ALL_TOOLS=1)
cmake .. -DCMAKE_INSTALL_PREFIX=..
make -j4
make test
make install
```


### Docker image

To build a docker image with overlap_checker and occ_faceter installed (using the Dockerfile in this directory), run:

```shell
docker build --pull --tag=occ_facet_geom .
```

Then to use the image, with data in the current directory appearing in /data:
```shell
docker run -v "$PWD:/data" -it occ_facet_geom
```

For development, use instead:

```shell
docker build --tag=inner --target inner .
```


### Additional make_watertight tests

There are tests under development which use [overlap_checker](https://github.com/ukaea/overlap_checker) tools to prepare test geometry, and then run the DAGMC make_watertight tool on the output of occ_faceter, with the expectation that make_watertight could highlight potential problems with occ_faceter's output.  See [prepare-test-geometries.sh](src/test/prepare-test-geometries.sh) and [test-with-make_watertight.sh](src/test/test-with-make_watertight.sh).


### Usage

After installation the bin directory should contain occ_faceter and other tools. To use occ_faceter you will need a CAD geometry,
either in BREP format (output from PPP, with an associated JSON file), or a JSON list of STEP files (with associated materials).

The example below assumes you have two files which were output by PPP, `test.brep` and the associated `test_metadata.json`.
The default values are used for faceting tolerance (relative to edge size) and output filename.  Run `./occ_faceter --help` to view the options and their default values.

mbconvert is used to convert the new h5m file into an vtk file which can then be visualized in Paraview (or Visit).  occ_faceter's optional add_mat_ids flag is used, so that different materials can be given different colours in Paraview.

```
cd ../bin
./occ_faceter --add_mat_ids test.brep
mbconvert dagmc_not_watertight.h5m test.vtk
paraview test.vtk
```
