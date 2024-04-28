# An Exact Algorithm for Disaster Resilient Planar Straight-Line Graph Augmentation
This repository contains a C++ implementation of an exact algorithm for solving a disaster resilient graph augmentation problem. 

## Dependencies
The solver depends upon [CGAL](doc.cgal.org) and the [Gurobi solver](https://www.gurobi.com/documentation/), as well as [boost](https://www.boost.org/doc/) (which is a dependency of CGAL).
Installation instructions can be found at the respective websites.
Gurobi version 11.0 is used, and the `makefile` will have to be modified to accommodate a different version. 
The makefile also assumes the definition of the `$GUROBI_HOME` shell variable. 

## Compilation

Compilation should work simply by running
```sh
make
```
assuming installation of the relevant dependencies. 
Installation has been tested on linux.
The directory can be cleaned of object files and the executable with `make clean`.

## Usage
Some example usage is given here, with more help available by running `./main --help`.

To run a batch solve on the instances listed in `instances/main-instances.txt`, running each instance with disasters incremented by 0.3, storing the results in the default `out.csv`, one would run
```sh
./main -i instances/main-instances.txt -d 0.3
```

To run a single instance with the disaster radius as specified in the input file, logging information to `log.txt`, one could run 
```sh
./main -s instance/network-deluauny-n-20-i-0.txt -l > log.txt
```