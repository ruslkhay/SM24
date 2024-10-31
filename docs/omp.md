# OpenMP implementation

The described problem requires calculation on plane. Therefor we are obliged to
iterate with double `for` loops.

## Usage in code

Pure OpenMP code is located in `src/openmp.cpp`. 

Calculation of $a_{ij}$, $b_{ij}$, $F_{ij}$, $w_{ij}$ requaries double 
for-loop. Therefor they were provided with OpenMP.

Used constructions:

1. **omp_set_dynamic(0)** \
This function is called to disable dynamic adjustment of the number of threads.
When you set it to 0, it means that the OpenMP runtime should not change the
number of threads during execution based on runtime conditions. By default, the
OpenMP runtime could adjust thread counts dynamically to optimize performance
based on workload, but disabling this can lead to more predictable behavior,
particularly beneficial for applications where workload is known in advance and
is stable.

2. **omp_set_num_threads(N)**  \
This function is used to set the number of threads to be used for parallel
sections of the code. Here, N would typically be an integer representing the
number of threads to use for parallel execution. By calling this
function, you can control parallelism manually rather than relying on OpenMP's
default behavior. 

Used pragmas:

1. **omp parallel for collapse(2)** \
This pragma is used to denote that the following loop (or loops) should be
executed in parallel.
The collapse(2) clause indicates that OpenMP should treat nested loops as 
a single loop with a larger number of iterations. This is useful
when you want to increase the granularity of the parallelism by combining the
iterations of multiple loops into a single set of iterations. For instance,
two loops that each iterate over M and N iterations respectively,
collapsing them will treat the combined iterations as M * N total iterations
for distributing among available threads.

2. **omp parallel for reduction(+ : res)** \
This pragma permits parallel execution of the following for loop, where the
results from each thread (or iteration) are combined at the end.
The reduction(+ : res) clause specifies that the variable res should be reduced
using the addition operator (+). Each thread computes its own partial sum, and
at the end of the parallel region, these partial sums are added together to
yield the final result stored in res. The reduction is crucial because it
handles potential data races on res by ensuring that each thread has its
private version of res that is then combined safely.

## Performance boost

Linear solution took 1803410 $\mu s$ (microseconds). 

| Threads   | Grid size  (M x N)  | Iter  | CPU Time ($\mu s$) | Boost %  |
|---|---|---|---|---|
| 1 | 40 x 40 | 100 000  | 2052314  | 87,9  |
| 4 |  40 x 40 | 100 000  | 1132357 | 159,3  |
| 16 |  40 x 40 | 100 000  | 2363495 | 76,3  |
