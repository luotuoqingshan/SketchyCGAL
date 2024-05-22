graphs = ["G$i" for i = 1:10]
seed = 0
R = 10
tol = 0.01
problem = "LovaszTheta"

open(homedir()*"/SketchyCGAL/test_$problem.txt", "w") do io
    for graph in graphs
        # make sure we warmup
        println(io, "ulimit -d $((16 * 1024 * 1024));"* 
        "matlab -singleCompThread -batch \"cd ~/SketchyCGAL/;"* 
        "Test_$(problem)_CGAL_PD('G1', $seed, $R, $tol);"*
        "Test_$(problem)_CGAL_PD('$graph', $seed, $R, $tol);\"")
    end
end