# OCC Faceter

OCC Faceter is a coolection of tools used for the production, modification and serialisation of geometries for [DAGMC](https://svalinn.github.io/DAGMC/)


### Installation

OCC Faceter requires [CGAL](https://cgal.org/), [MOAB](https://press3.mcs.anl.gov/sigma/moab-library/) and [OCE](https://github.com/tpaviot/oce) or [OCC](https://www.opencascade.com) to run

Install the dependencies and start using the tool.

```sh
$ sudo apt-get install libcgal-dev
$ sudo add-apt-repository ppa:freecad-maintainers/freecad-stable
$ sudo apt-get install libocc*-dev
$ sudo apt-get install occ*
$ sudo apt-get install libtbb-dev
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_INSTALL_PREFIX=..
$ make
$ make install
```

### Usage

After installation the $bin$ directory should contain $steps2h5m$ and other tools. To use $steps2h5m$ you will need a json file containing a list of step files and their associated material names. Below is an example JSON_FILE.

$
[
    {
        "material": "steel",
        "filename": "tube.stp"
    },
    {
        "material": "water",
        "filename": "coolant.stp"
    }
]
$

The general useage of steps2h5m requires three arguments.

```
$ ./steps2h5m JSON_FILE TOLERANCE H5M_FILE
```

Assuming the JSON_FILE is called geometry_desciption.json and the tolerance is 0.1 and the output file is dagmc.h5m then the command would be


```
$ ./steps2h5m geometry_desciption.json 0.1 dagmc.h5m
```





