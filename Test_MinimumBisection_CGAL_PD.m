function out = Test_MinimumBisection_CGAL_PD(varargin)
    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'seed', 0, @isnumeric);
    addOptional(p, 'R', 10, @isnumeric);
    addOptional(p, 'tol', 0.01, @isnumeric)

    parse(p, varargin{:});

    seed = p.Results.seed; % random seed
    minbisec_data = p.Results.graph;
    R = p.Results.R; % rank/sketch size parameter
    tol = p.Results.tol;

    % maxcut_data = 'DIMACS10/belgium_osm';
    
    fprintf("Solving Minimum Bisection SDP for %s\n", minbisec_data);
    %% Preamble
    rng(seed,'twister');
    %% Please update these paths before running 
    addpath /homes/huan1754/SketchyCGAL/utils;
    addpath /homes/huan1754/SketchyCGAL/solver;
    
    %% Load data
    
    %% Modify the path before running
    data = load(['~/datasets/graphs/MinimumBisection/', minbisec_data, '.mat']);
    A = data.A;
    
    n = size(A,1);
    C = spdiags(A*ones(n,1),0,n,n) - A;
    C = 0.5*(C+C'); % symmetrize if not symmetric
    C = 0.25.*C;
    e = ones(n, 1);
    
    clearvars Problem;
    clearvars data;
    
    %% Construct the Black Box Oracles for MaxCut SDP
    
    Primitive1 = @(x) C*x;
    Primitive2 = @(y,x) (y(1:end-1, 1).*x +y(end, 1)*(e'*x) * e);
    Primitive3 = @(x) [x.^2; (e'*x).^2];

    a = n; % trace bound
    b = [ones(n,1); 0.0];
    
    % Compute scaling factors
    SCALE_X = 1/n;
    SCALE_C = 1/norm(C,'fro');
    
    %% Solve using SketchyCGAL
    
    beta0 = 1; % we didn't tune - choose 1 - you can tune this!
    K = inf;
    maxit = 1e6; % limit on number of iterations
    
    timer = tic;
    cputimeBegin = cputime;

    [out, U, D, y, AX, pobj] = CGAL( n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
        'FLAG_MULTRANK_P1',true,... % This flag informs that Primitive1 can be applied to find AUU' for any size U. 
        'FLAG_MULTRANK_P3',true,... % This flag informs that Primitive1 can be applied to find (A'y)U for any size U.
        'SCALE_X',SCALE_X,... % SCALE_X prescales the primal variable X of the problem
        'SCALE_C',SCALE_C,... % SCALE_C prescales the cost matrix C of the problem
        'stoptol', tol,...
        'evalsurrogategap', true); % algorithm stops when 1e-3 accuracy is achieved
                         
    cputimeEnd = cputime;
    totalTime = toc(timer);
    
    out.totalTime = totalTime;
    out.totalCpuTime = cputimeEnd - cputimeBegin;
    
    %% Evaluate errors
    
    
    cutvalue = Inf;
    for repeat = 1:100
        x = U*randn(size(U, 2), 1);
        [~, ord] = sort(x);
        partition_vec = zeros(n, 1);
        for i = 1:n
            if i <= n/2
                partition_vec(ord(i)) = 1;
            else
                partition_vec(ord(i)) = -1;
            end
        end
        rankvalue = (partition_vec'*(C*partition_vec));
        cutvalue = min(cutvalue, rankvalue);
    end
    
    out.cutvalue = cutvalue;
    out.primalObj = pobj;
    out.primalFeas = norm(AX-b)/(1+norm(b));
    disp(out.primalObj);
    disp(out.cutvalue);
    
    %% Save results
    
    
    if ~exist(['~/SDPLR.jl/output/MinimumBisection/',minbisec_data,...
         '/SketchyCGAL'],'dir') 
        mkdir(['~/SDPLR.jl/output/MinimumBisection/',minbisec_data,...
         '/SketchyCGAL']); 
    end
    save(['~/SDPLR.jl/output/MinimumBisection/',minbisec_data,...
        '/SketchyCGAL/SketchyCGAL-R-', num2str(R), '-seed-',...
        num2str(seed), '-tol-', num2str(tol), '.mat'],'out','-v7.3');
end
%% Last edit: Alp Yurtsever - July 24, 2020