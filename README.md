# Introduction  

This code repo provides implementation for applying 
[SketchyCGAL](https://github.com/alpyurtsever/SketchyCGAL) to 
three more SDPs, Minimum Bisection, Lovasz Theta and Cut Norm.  It also provides some easy-to-use batch test tool for massive benchmarking.

## Data
Data used in our experiments can be downloaded from [University of Florida Sparse Matrix Collection](https://sparse.tamu.edu/) or [SNAP](https://snap.stanford.edu/data/index.html). We provide one example script for downloading DIMACS10 data, see `download.sh`.
Here we also include some preprocessed toy Gset graphs for trial usage. Basically each 
mat file contains one adjacency matrix.

In general, data loading can be easily customized, just modify the data loading part in each 
code file (e.g. line 29-30 of Test_LovaszTheta_CGAL_PD.m).  

## Running  
For a specific SDP problem, its code is in the corresponding Test_*Problem*_CGAL_PD.m file. You can either call the function directly inside 
matlab or run the batch test script from the shell.

To call the function directly, make sure you are in the `SketchyCGAL/` folder, 
and feed the path to this folder into the function, e.g. 
```
cd ~/SketchyCGAL/
Test_CutNorm_CGAL_PD(path2folder='~/')
```
You may need to modify the path a bit to fit your case.

For batch testing, you can run `julia gen_test.jl` to generate a txt file 
where each line corresponds to one problem instance (install julia if needed). You may want to 
change the parameters on the top of `gen_test.jl` to fit your need. By default we allow solver 16GB RAM access.
For fair of comparison, we recommend you to **disable multithreading**. For example, you can add
the following three lines into your `.bashrc`, and also call matlab with flag 
`-singleCompThread` (which is automatically set in `gen_test.jl`). 
```
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```
After you run `gen_test.jl` and the batch test txt file is generated, 
you can parallel the benchmarking via using GNU parallel, e.g. 
```
cat test_MaxCut.txt | parallel --jobs 10 --timeout 28800 {}
```
Here `--timeout` specifies the time limit, and `--jobs` means the number of 
parallel jobs you want to run. Keep in mind that you should not set `--jobs` 
larger than the number of cores you have, which will downgrade the performance.

# Acknowledgement
Thanks for the great toolbox `SketchyCGAL` which provides an easy-to-use interface and allows quick development of these experiments. 