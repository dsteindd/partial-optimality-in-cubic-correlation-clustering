# Code for the Paper "Partial Optimality in Cubic Correlation Clustering"

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23643382.svg)](https://doi.org/10.5281/zenodo.23643382)

This is the code for [Partial Optimality in Cubic Correlation Clustering](https://proceedings.mlr.press/v202/stein23a.html), published at ICML 2023 (also on [arXiv](https://arxiv.org/abs/2302.04694)).

## Requirements
- CMake >= 3.10 and a C++17 compiler
- the BOOST Graph Library (BGL), see https://www.boost.org/doc/libs/1_81_0/libs/graph/doc/index.html
- Python 3 with pandas, to aggregate the results (see [Aggregation](#aggregation))

The random instances are drawn with the distributions of the C++ standard library, which differ between implementations.
The results of the paper are reproduced exactly with GCC and its standard library libstdc++ (tested with GCC 16.2).
Since GCC 16, `std::generate_canonical` follows P0952R2, which changes the random numbers.
The `CMakeLists.txt` therefore defines `_GLIBCXX_USE_OLD_GENERATE_CANONICAL`, which restores the previous definition.
With other standard libraries, e.g. libc++ of Clang, the instances and hence the numbers differ.

## How to build?

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make equilateral-triangles-persistency
make predefined-partition-persistency
make equilateral-triangles-instance-sizes
make predefined-partition-instance-sizes

./equilateral-triangles-persistency
./predefined-partition-persistency
```

## How to run?
Depending on what problem instances you want to run, execute one of the following commands.
Note, that we run the computation 30 times with different random seeds for each experiment.
The drivers append to existing output files, so delete the output directory before running a driver again.

### Equilateral Triangle Instances / Geometric Instances

#### Joint Application
The default is to apply all our persistency conditions jointly to the instance of equilateral triangles with 45 nodes.
For this default, run:
```shell
./equilateral-triangles-persistency -a
```
or 
```shell
./equilateral-triangles-persistency --all
```
Alternatively, specify the number of elements of the geometric instance in the following way. Run
```shell
./equilateral-triangles-persistency -a <n>
```
where `<n>` should be replaced by a natural number greater or equal to 3. 
Here, `<n>` is a design parameter and leads to a geometric instance with a total of `N= 9n - 9` elements, i.e. by default `<n>=6`.


#### Separate Application
The default is to apply all our persistency conditions separately to the instance of equilateral triangles with 45 nodes.
Run:
```shell
./equilateral-triangles-persistency -i
```
or
```shell
./equilateral-triangles-persistency --individual
```
Alternatively, specify the number of elements of the geometric instance in the following way. Run
```shell
./equilateral-triangles-persistency -i <n>
```
where `<n>` should be replaced by a natural number greater or equal to 3.
Here, `<n>` is a design parameter and leads to a geometric instance with a total of `N= 9n - 9` elements, i.e. by default `<n>=6`.


### Partition Instances
#### Joint Application
The default is to apply all our persistency conditions jointly to the partition instances with 48 nodes.
```shell
./predefined-partition-persistency -a
```
or 
```shell
./predefined-partition-persistency --all
```
To specify the size of the instance run
```shell
./predefined-partition-persistency -a <n>
```
where `<n>` is a natural number and leads to an instance of the partition dataset with `N = 8 <n>` elements.

#### Separate Application
The default is to apply all our persistency conditions separately to the partition instances with 48 nodes.
```shell
./predefined-partition-persistency -i
```
or
```shell
./predefined-partition-persistency --individual
```
To specify the size of the instance run
```shell
./predefined-partition-persistency -i <n>
```
where `<n>` is a natural number and leads to an instance of the partition dataset with `N = 8 <n>` elements.


### Instance Sizes
To apply all our persistency conditions jointly to instances of growing size, run
```shell
./equilateral-triangles-instance-sizes
./predefined-partition-instance-sizes
```
The geometric instances are run for `sigma` in `{0.01, 0.06, 0.1}` with `N = 9n - 9` elements for `n = 2, ..., 13`.
The partition instances are run for `alpha` in `{0.0, 0.4, 0.5, 0.65}` and `beta = 0.5` with `N = 8n` elements for `n = 1, ..., 13`.
By default, each configuration is run with 30 random seeds.
Use `--seeds <number>` to change the number of seeds and `--max-n <n>` to stop at a smaller instance size, e.g.
```shell
./predefined-partition-instance-sizes --seeds 10 --max-n 6
```


## Other Configurations
If you want to try different configurations, you can change the parameters in the files in 

- `./src/equilateral-triangles/persistency.cxx`
- `./src/predefined-partition/persistency.cxx`
- `./src/equilateral-triangles/instance-sizes.cxx`
- `./src/predefined-partition/instance-sizes.cxx`

## Output
### Joint Application
The output of the joint application of our partial optimality criteria is a directory named

- `partition_<numNodes>nodes_30seeds` for the partition instances
- `equilateral_<numNodes>nodes_30seeds_30distanceOfCostZero_4sigma` for the geometric instances

These contain the number of eliminated nodes/variables after each step of the joint application of our criteria.


### Separate Application
The output of the separate application of our partial optimality criteria is a directory named

- `partition_individual_<numNodes>nodes_30seeds` for the partition instances
- `equilateral_individual_<numNodes>nodes_30seeds_30distanceOfCostZero_4sigma` for the geometric instances

These contain, for each partial optimality criterion, a `.csv` file containing the number eliminated nodes/variables.

### Instance Sizes
The output of the instance size experiments is a directory named

- `partition_instance_sizes_<numSeeds>seeds` for the partition instances
- `equilateral_instance_sizes_<numSeeds>seeds_30distanceOfCostZero_4sigma` for the geometric instances

These contain one `.csv` file per instance size with the number of eliminated nodes/variables after each step of the joint application of our criteria.

## Aggregation
The script `scripts/aggregate.py` turns the output directories into the data files of the figures in the paper, e.g.
```shell
python3 ../scripts/aggregate.py equilateral_45nodes_30seeds_30distanceOfCostZero_4sigma partition_48nodes_30seeds --out aggregated
```
It recognizes the experiment from the name of each directory and writes, for every configuration, statistics over the seeds:

- the percentage of fixed variables, `100 * (initial - final) / initial`, where `final` is the number of remaining variables after the last step (for `tripletCut` also the percentage of fixed triples),
- the mean, standard deviation, median, lower and upper quartile, minimum and maximum of these percentages,
- the same statistics of the runtime in milliseconds, i.e. the sum of the `runTime` column of a run.

## Citation
If you use this code, please cite the paper:
```bibtex
@inproceedings{stein-2023-partial,
    author = {David Stein and Silvia Di Gregorio and Bjoern Andres},
    title = {Partial Optimality in Cubic Correlation Clustering},
    booktitle = {ICML},
    year = {2023},
    url = {https://proceedings.mlr.press/v202/stein23a.html},
}
```
The code is archived on Zenodo under the DOI [10.5281/zenodo.23643382](https://doi.org/10.5281/zenodo.23643382), which always refers to the latest version.

## License
This code is licensed under the MIT License, see [LICENSE](LICENSE).
The exception is `include/andres/partition.hxx`, which is taken from the [graph library of Bjoern Andres](http://www.andres.sc/graph.html) and distributed under the license stated in that file.
