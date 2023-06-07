# Code for the Paper "Partial Optimality in Cubic Correlation Clustering"

This is the code for [Partial Optimality in Cubic Correlation Clustering](https://arxiv.org/abs/2302.04694), accepted at ICML'23.

## Requirements
You need to have the BOOST Graph Library (BGL) installed. See here: https://www.boost.org/doc/libs/1_81_0/libs/graph/doc/index.html

## How to build?

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make equilateral-triangles-persistency
make predefined-partition-persistency

./equilateral-triangles-persistency
./predefined-partition-persistency
```

## How to run?
Depending on what problem instances you want to run, execute one of the following commands.
Note, that we run the computation 30 times with different random seeds for each experiment.

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


## Other Configurations
If you want to try different configurations, you can change the parameters in the files in 

- `./src/equilateral-triangles/persistency.cxx`
- `./src/predefined-partition/persistency.cxx`

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
