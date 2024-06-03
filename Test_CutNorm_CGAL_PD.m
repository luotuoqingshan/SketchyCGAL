function out = Test_CutNorm_CGAL_PD(varargin)
    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'seed', 0, @isnumeric);
    addOptional(p, 'R', 10, @isnumeric);
    addOptional(p, 'tol', 0.01, @isnumeric);
    addOptional(p, 'path2folder', '~/', @ischar)

    parse(p, varargin{:});

    seed = p.Results.seed; % random seed
    cutnorm_data = p.Results.graph;
    R = p.Results.R; % rank/sketch size parameter
    tol = p.Results.tol; % stopping tolerance
    path2folder = p.Results.path2folder; % path to the folder of this repo


    fprintf("Solving CutNorm SDP for %s\n", cutnorm_data);
    %% Preamble
    rng(seed,'twister');

    %% Please update these paths before running 
    addpath utils;
    addpath solver;

    %% Load data

    data = load([path2folder, 'SketchyCGAL/data/CutNorm/', cutnorm_data, '.mat']);
    A = data.A;
    [m, n] = size(A);

    C = [sparse(m, m) A; A' sparse(n, n)];
    C = 0.5*C; % symmetrize if not symmetric
    C = -C; % negate to convert to minimization problem
    % Interestingly, since cut norm is a norm
    % cut norm of A and -A are the same
    % so this negation is not necessary

    clearvars Problem;
    clearvars data;

    %% Construct the Black Box Oracles for CutNorm SDP
    Primitive1 = @(x) C*x;
    Primitive2 = @(y,x) y.*x;
    Primitive3 = @(x) sum(x.^2,2);
    a = m+n; % trace bound
    b = ones(m+n,1);

    % Compute scaling factors
    SCALE_X = 1/(m+n);
    SCALE_C = 1/norm(C,'fro');

    %% Solve using SketchyCGAL

    beta0 = 1; % we didn't tune - choose 1 - you can tune this!
    K = inf;
    maxit = 1e6; % limit on number of iterations

    timer = tic;
    cputimeBegin = cputime;

    [out, ~, ~, ~, AX, pobj] = CGAL(m+n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
        'FLAG_MULTRANK_P1',true,... % This flag informs that Primitive1 can be applied to find AUU' for any size U. 
        'FLAG_MULTRANK_P3',true,... % This flag informs that Primitive1 can be applied to find (A'y)U for any size U.
        'SCALE_X',SCALE_X,... % SCALE_X prescales the primal variable X of the problem
        'SCALE_C',SCALE_C,... % SCALE_C prescales the cost matrix C of the problem
        'stoptol',tol,...
        'evalsurrogategap', true); 

    cputimeEnd = cputime;
    totalTime = toc(timer);

    out.totalTime = totalTime;
    out.totalCpuTime = cputimeEnd - cputimeBegin;

    out.primalObj = pobj;
    out.primalFeas = norm(AX-b)/(1+norm(b));

    %% Save results
    if ~exist([path2folder, 'SketchyCGAL/output/CutNorm/',cutnorm_data],'dir') 
        mkdir([path2folder, 'SketchyCGAL/output/CutNorm/',cutnorm_data]); 
    end
    save([path2folder, 'SketchyCGAL/output/CutNorm/',cutnorm_data,...
        '/SketchyCGAL-R-', num2str(R), '-seed-',...
        num2str(seed), '-tol-', num2str(tol), '.mat'],'out','-v7.3');
end