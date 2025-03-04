# mpi

This repo contains several mpi apps and algorithm based fault tolerant apps (mm, lu, kmeans).

- CRC32 is used to verify the integrity of communication data.
- Hamming code (SECDED) is used to correct one-bit errors and detect two-bit errors in data blocks. 

## MM
- Compile
```
make mm
make mm_abft
make mm_comp
```
- Run
```
mpirun -np <num_of_procs> bin/mm <matrix_a_filename> <matrix_b_filename>
mpirun -np <num_of_procs> bin/mm_abft <matrix_a_filename> <matrix_b_filename>
mpirun -np <num_of_procs> bin/mm_comp <matrix_a_filename> <matrix_b_filename> [CT]
```
- A small program [gen_matrices.py](src/mm/gen_matrices.py) helps to generate a matrix. 
- Requirement
    - m_1->cols == m_2->rows
    - m_1->rows % num_worker == 0

## LU
- Compile
```
make lu
make lu_abft
make lu_comp
```
- Run
```
mpirun -np <num_of_procs> bin/lu <matrix_size>
mpirun -np <num_of_procs> bin/lu_abft <matrix_size>
mpirun -np <num_of_procs> bin/lu_comp <matrix_size> [CT]
```

## K-means
- Compile
```
make kmeans
make kmeans_abft
make kmeans_comp
```
- Run
```
mpirun -np <num_of_procs> bin/kmeans [clusters] [max_iterations] [datafile] # (default) 100 1000 ./data/obs_info.txt
mpirun -np <num_of_procs> bin/kmeans_abft [clusters] [max_iterations] [datafile] # (default) 100 1000 ./data/obs_info.txt
mpirun -np <num_of_procs> bin/kmeans_comp [CT] # CT, clusters, max_iterations, datafile are set in include/dataCompression.h
```

## Himeno
- Parameter Setting
```
cd src/himeno
./paramset.sh M 1 1 2 # generate new param.h
```
- Compile
```
make himeno
make himeno_comp
```
- Run
```
mpirun -np <num_of_procs> bin/himeno
```

## Graph500
- Compile
```
cd src/graph500/mpi/
make
```
- Run
```
cd src/graph500/mpi/
./graph500_mpi_simple <scale> [edge-factor-default-16]
```

## PingPong
- Compile
```
make pingpong_float_comp
make pingpong_double_comp
```
- Run
```
mpirun -np 2 ./bin/pingpong_float_comp
mpirun -np 2 ./bin/pingpong_double_comp
```