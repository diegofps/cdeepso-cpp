# About this project

This is a fast c++ implementation of the C-DEEPSO algorithm. It was inspired by the [Matlab implementation](https://www.mathworks.com/matlabcentral/fileexchange/66980-c-deepso-canonical-differential-evolutionary-particle-swar) and the algorithm's [original paper](https://link.springer.com/article/10.1007/s10489-018-1167-5).

# Dependencies

1. wup - a general purpose library with many of my reusable utilities.

2. clang++ - used to compile the program.

# Getting started

Clone the repository and its submodule (wup).

```shell
git clone --recurse-submodules git@github.com:diegofps/cdeepso-cpp.git
```

Access the src folder and type make to build the example.

```shell
cd cdeepso-cpp/src
make
```

Running it.

```shell
# Default parameters
./main

# rastrigin function
./main -eval ras

# rosenbrock function
./main -eval ros

# griewank function
./main -eval gri

# Specify the number of threads, e.g., 16
./main -threads 16

# Same number of threads as CPU cores
./main -threads 0

# Max fitness evals
./main -maxFitEval 100000

# Max generations
./main -maxGen 50000

# Population size
./main -popSize 50

# Dimensinality
./main -dims 50

# Mutation rate
./main -mutationRate 0.5

# Communication probability
./main -communicationProbability 2.0

# Max particle velocity
./main -maxVelocity 2.0

# The number of times CDEEPSO will run (to obtain mean and average results)
./main -maxRun 50
```

# Performance results

# Results 2 (2023)

The following results were obtained in a machine running Ubuntu 22.04.2 LTS, AMD Threadripper 2990WX and 128GB of RAM. The parameters used are shown in the command lines bellow the table.

| Function | Mean fitness (std) | Execution time in ms (std) |
| -------- | ------------------ | -------------------------- |
| ras      |  4.94 (4.59)       | 498.23 (20.72)             |
| ros      | 95.00 (1.27)       | 312.46 (20.74)             |
| gri      |  0.00 (0.00)       | 454.37 (20.56)             |

```shell
# RAS
./main \
    -maxFitEval 100000 \
    -maxGen 50000 \
    -popSize 50 \
    -dims 100 \
    -memGBestSize 5 \
    -eval ras \
    -memStrategy POS_MEM \
    -deType RAND \
    -mutationRate 0.5 \
    -communicationProbability 0.1 \
    -printConvergenceResults 0 \
    -maxRun 1000 \
    -threads 64

# ROS
./main \
    -maxFitEval 100000 \
    -maxGen 50000 \
    -popSize 50 \
    -dims 100 \
    -memGBestSize 5 \
    -eval ros \
    -memStrategy POS_MEM \
    -deType RAND \
    -mutationRate 0.5 \
    -communicationProbability 0.1 \
    -printConvergenceResults 0 \
    -maxRun 1000 \
    -threads 64

# GRI
./main \
    -maxFitEval 100000 \
    -maxGen 50000 \
    -popSize 50 \
    -dims 100 \
    -memGBestSize 5 \
    -eval gri \
    -memStrategy POS_MEM \
    -deType RAND \
    -mutationRate 0.5 \
    -communicationProbability 0.1 \
    -printConvergenceResults 0 \
    -maxRun 1000 \
    -threads 64
```

# Results 1 (2020)

The following results were obtained in a machine running Ubuntu 20.04.4 LTS, AMD Threadripper 2990WX and 128GB of RAM. The number of threads was 64 (-threads), the particles had 50 dimensions (-dims), and we executed CDEEPSO 100 times (-maxRun).

| Function | Mean fitness | Std fitness | Execution time (ms) |
| -------- | ------------ | ----------- | ------------------- |
| ras      | 0            | 0           | 470.14              |
| ros      | 43.46        | 2.31        | 317.91              |
| gri      | 0            | 0           | 448.81              |
