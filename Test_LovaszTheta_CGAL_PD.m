function out = Test_LovaszTheta_CGAL_PD(varargin)
    p = inputParser;
    addOptional(p, 'graph', 'G1', @ischar);
    addOptional(p, 'seed', 0, @isnumeric);
    addOptional(p, 'R', 10, @isnumeric);
    addOptional(p, 'tol', 0.01, @isnumeric);

    parse(p, varargin{:});

    seed = p.Results.seed; % random seed
    lovasztheta_data = p.Results.graph;
    R = p.Results.R; % rank/sketch size parameter
    tol = p.Results.tol; % stopping tolerance


    fprintf("Solving Lovasz Theta SDP for %s\n", lovasztheta_data);
    %% Preamble
    rng(seed,'twister');
    %% Please update these paths before running 
    addpath /homes/huan1754/SketchyCGAL/utils;
    addpath /homes/huan1754/SketchyCGAL/solver;

    %% Load data

    %% Modify the path before running
    % We use the same graph for Max Cut and Lovasz Theta
    data = load(['~/datasets/graphs/MaxCut/', lovasztheta_data, '.mat']);
    A = data.A;

    n = size(A,1);
    e = ones(n,1);

    clearvars Problem;
    clearvars data;

    %% Construct the Black Box Oracles for MaxCut SDP

    function out = Astary_x(y, x, i, j)
        n = length(x);
        out = zeros(n, 1); 
        for k = 1:length(i)
            if i(k) == j(k)
                out(i(k)) = out(i(k)) + y(k) * x(i(k));
            else
                out(i(k)) = out(i(k)) + y(k) * x(j(k)) / 2;
                out(j(k)) = out(j(k)) + y(k) * x(i(k)) / 2;
            end
        end
        e = ones(n, 1);
        out = out + y(end) .* x;
    end

    function out = A_xtx(x, i, j)
        n = length(x);
        m = length(i);
        out = zeros(m+1, 1);
        for k = 1:length(i)
            out(k) = x(i(k)) * x(j(k));
        end
        out(end) = sum(sum(x.^2));
    end

    [i, j, ~] = find(triu(A));
    m = length(i);
    Primitive1 = @(x) -sum(e.*x).*e;
    Primitive2 = @(y,x) Astary_x(y, x, i, j);
    Primitive3 = @(x) A_xtx(x, i, j);
    a = 1; % trace bound
    b = zeros(m+1,1);
    b(end) = 1;

    % Compute scaling factors
    SCALE_X = 1;
    SCALE_C = 1/n;

    %% Solve using SketchyCGAL

    beta0 = 1; % we didn't tune - choose 1 - you can tune this!
    K = inf;
    maxit = 1e6; % limit on number of iterations

    timer = tic;
    cputimeBegin = cputime;

    [out, U, D, y, AX, pobj] = CGAL(n, Primitive1, Primitive2, Primitive3, a, b, R, maxit, beta0, K, ...
        'FLAG_MULTRANK_P1',true,... % This flag informs that Primitive1 can be applied to find AUU' for any size U. 
        'FLAG_MULTRANK_P3',true,... % This flag informs that Primitive1 can be applied to find (A'y)U for any size U.
        'SCALE_X',SCALE_X,... % SCALE_X prescales the primal variable X of the problem
        'SCALE_C',SCALE_C,... % SCALE_C prescales the cost matrix C of the problem
        'stoptol',tol,...
        'evalsurrogategap', true,...
        'carefulstopping', true); 
        % SketchyCGAL stops quickly on lovasz theta problems 
        % so we use careful stopping to ensure we get an accurate duality bound 

    cputimeEnd = cputime;
    totalTime = toc(timer);

    out.totalTime = totalTime;
    out.totalCpuTime = cputimeEnd - cputimeBegin;
    out.primalObj = pobj;
    out.primalFeas = norm(AX-b)/(1+norm(b));

    %% Save results

    if ~exist(['~/SketchyCGAL/output/LovaszTheta/',lovasztheta_data],'dir') 
        mkdir(['~/SketchyCGAL/output/LovaszTheta/',lovasztheta_data]); 
    end
    save(['~/SketchyCGAL/output/LovaszTheta/',lovasztheta_data,...
        '/SketchyCGAL-R-', num2str(R), '-seed-',...
        num2str(seed), '-tol-', num2str(tol), '.mat'],'out','-v7.3');
end
%% Last edit: Alp Yurtsever - July 24, 2020