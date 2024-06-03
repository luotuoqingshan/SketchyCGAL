# modify following keywords for your need
graphs = ["G$i" for i = 1:10] # list of graphs to benchmark
seed = 0 # random seed
R = 10 # sketch size
tol = 0.01 # tolerance for primal infeasiblity and suboptimality(relative)
problem = "MinimumBisection" # problem you want to solve
path2folder = "~/" # path to the folder SketchyCGAL/

##
open(homedir()*"/SketchyCGAL/test_$problem.txt", "w") do io
    for graph in graphs
        # make sure we warmup
        println(io, "ulimit -d $((16 * 1024 * 1024));"* 
        "matlab -singleCompThread -batch \"cd ~/SketchyCGAL/;"* 
        "Test_$(problem)_CGAL_PD('G1', $seed, $R, $tol, '$path2folder');"*
        "Test_$(problem)_CGAL_PD('$graph', $seed, $R, $tol, '$path2folder');\"")
    end
end