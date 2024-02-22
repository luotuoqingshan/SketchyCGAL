%% Test setup for Primal Dual Convergence (MaxCut SDP) - Solved with SketchyCGAL

%% Choose data
% NOTE: You need to download data from GSET and locate them to under the
% "FilesMaxCut/data/G/" folder (resp. DIMACS10, "FilesMaxCut/data/DIMACS10/").

for gid = 1:67
    maxcut_data = ['G' num2str(gid)]; % you can choose other data files as well
    % maxcut_data = 'DIMACS10/belgium_osm';
    
    fprintf("Solving Minimum Bisection SDP for %s\n", maxcut_data);
    %% Preamble
    rng(0,'twister');
    addpath utils;
    addpath solver;
    
    %% Load data
    
    %load(['./FilesMaxCut/data/',maxcut_data]);
    A = readSMAT(['~/Gset/', maxcut_data, '.smat']);
    
    n = size(A,1);
    disp(size(A))
    C = spdiags(A*ones(n,1),0,n,n) - A;
    C = 0.5*(C+C'); % symmetrize if not symmetric
    C = 0.25.*C;
    d = sum(A, 2);
    
    clearvars Problem;
    
    %% Construct the Black Box Oracles for MaxCut SDP
    
    Primitive1 = @(x) C*x;
    Primitive2 = @(y,x) (y(1:end-1, 1).*x +y(end, 1)*(d'*x) * d);
    Primitive3 = @(x) [x.^2; (d'*x).^2];
    a = [0,1.05*n];
    b = [ones(n,1); 0.0];
    
    % Compute scaling factors
    SCALE_X = 1/n;
    SCALE_C = 1/norm(C,'fro');
    
    %% Solve using SketchyCGAL
    
    R = 10; % rank/sketch size parameter
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
        'stoptol',1e-1,...
        'evalsurrogategap', true); % algorithm stops when 1e-3 accuracy is achieved
                         
    cputimeEnd = cputime;
    totalTime = toc(timer);
    
    out.totalTime = totalTime;
    out.totalCpuTime = cputimeEnd - cputimeBegin;

    disp(out.totalTime);
    disp(out.totalCpuTime);
    
    %% Evaluate errors
    
    dobj = b'*y;
    n = size(C, 1);

    eigPrimitive = @(x) (C* x + y(1:end-1, 1).*x +y(end, 1)*(d'*x) * d); 

    disp(eigs(eigPrimitive, n, 1, 'smallestreal'));

    out.dimacs.err1 = norm(AX-b)/(1+norm(b));
    out.dimacs.err2 = 0; % this is theoretically 0 for CGAL
    out.dimacs.err3 = 0; % this is theoretically 0 for CGAL
    out.dimacs.err4 = max(-min(eigs(eigPrimitive, n, 1, 'smallestreal')),0)/(1+norm(C,'fro'));
    out.dimacs.err5 = (pobj+dobj)/(1+abs(pobj)+abs(dobj));
    out.dimacs.err6 = (pobj + y'*AX)/(1+abs(pobj)+abs(dobj));
    
    cutvalue = 0;
    for t = 1:R
        sign_evec = sign(U(:,t));
        rankvalue = -(sign_evec'*(C*sign_evec));
        cutvalue = max(cutvalue, rankvalue);
    end
    
    out.cutvalue = cutvalue;
    out.primalObj = pobj;
    out.primalFeas = norm(AX-b)/(1+norm(b));
    
    %% Save results
    
    if ~exist(['results/DualMinimumBisection/',maxcut_data],'dir'), mkdir(['results/DualMinimumBisection/',maxcut_data]); end
    save(['results/DualMinimumBisection/',maxcut_data,'/SketchyCGAL.mat'],'out','-v7.3');
    
end
%% Last edit: Alp Yurtsever - July 24, 2020