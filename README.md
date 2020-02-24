# OCC Faceter

OCC Faceter is a coolection of tools used for the production, modification and serialisation of geometries for [DAGMC](https://svalinn.github.io/DAGMC/)


### Installation

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
sudo apt install software-properties-common # for add-apt-repository
sudo add-apt-repository ppa:freecad-maintainers/freecad-stable # for occ
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
cmake .. -DCMAKE_INSTALL_PREFIX=..
make -j4
make test
make install
```

### Usage

After installation the bin directory should contain occ_faceter and other tools. To use occ_faceter you will need a CAD geometry saved in Step file format. Th example below assumes you have a Step file called test.step. The faceting tolerance is set as 1.e-3 and the output file is called test.h5m. As an optional step mbconvert is used to convert the new h5m file into an stl file which can be visualized in paraview or visit  

```
$ cd ../bin
$ ./occ_faceter test.step 1.e-3 test.h5m
$ mbconvert test.h5m test.stl
```



