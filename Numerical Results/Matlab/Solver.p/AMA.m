%% ---------------------------------------------------------
%% This is the implementation of AMA 
%% Reference:
%% Eric C. Chi & Kenneth Lange, Splitting methods for convex clustering, JCGS 2015
%% parameter:
%% lambda: hyper-parameter of the convex clustering model
%% mu: hyper-parameter for the solver
%% Algorithm Copyright (c) 2015 Eric C. Chi and Kenneth Lange
%% Code Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% ------------------------------------------------------------
function [Z, X, info] = AMA(Ainput,b,dim,weightVec,options,Z0)
   maxiter = 2000;  
   mu = 2/(dim.n+1);
   stoptol = 1e-6;
   rng('default');
   if isfield(options,'stoptol'); stoptol = options.stoptol; end
   if isfield(options,'maxiter');  maxiter = options.maxiter; end
   %% Initialization
   fprintf('\n *******************************************************');
   fprintf('******************************************');
   fprintf('\n \t\t   Fast AMA  for solving Clustering with  mu = %6.3f', mu);
   fprintf('\n ******************************************************');
   fprintf('*******************************************\n');
   fprintf(' problem size: d = %3.0f, n = %3.0f, E = %3.0f',dim.d, dim.n, dim.E);
   fprintf('\n ---------------------------------------------------');
   fprintf('\n  iter| primobj | dualobj | relative gap | Time(s)');
   Fnorm = @(x) mexFnorm(x); 
   t_start = clock();
   Amap = Ainput.Amap;
   ATmap = Ainput.ATmap;
   ATAmap = Ainput.ATAmap;
   Amapbvec = Amap(b);
   if (~exist('Z0','var'))
       Z = sparse(dim.d, dim.E);
   else
       Z = Z0;
   end
   X = ATmap(Z) + b;
   info.iter = 0;
   info.dualobj = 0;
   info.primobj = 0.5*Fnorm(X - b)^2 + weightVec*(sqrt(sum(Amap(X).^2)))';
   info.dualgap = abs(info.primobj - info.dualobj)/(1.0 + abs(info.primobj) + abs(info.dualobj));
   alpha_old = 1;
   fprintf('\n %5.0d| %- 5.4e | %- 5.4e | %- 5.4e | %- 5.4e ', info.iter, info.primobj, info.dualobj, info.dualgap, etime(clock(), t_start));
   while (info.dualgap > stoptol)
       info.iter = info.iter + 1;
       Delta = ATmap(Z);
       gvec = Amapbvec + Amap(Delta);
       proj_input = Z - mu*gvec;
       Z = proj_l2(proj_input, weightVec);
       %% compute duality gap
       ATZ = ATmap(Z);
       ATZ_norm = Fnorm(ATZ);
       info.dualobj = -0.5*ATZ_norm^2 - sum(sum(ATZ.*b));
       %% Recovery primal variable
       X = b + ATZ;
       info.primobj = 0.5*ATZ_norm^2 + weightVec*(sqrt(sum(Amap(X).^2)))';
       info.dualgap = abs(info.primobj - info.dualobj)/(1.0 + abs(info.primobj) + abs(info.dualobj));
       fprintf('\n %5.0d| %- 5.4e | %- 5.4e | %- 5.4e | %- 5.4e ', info.iter, info.primobj, info.dualobj, info.dualgap, etime(clock(), t_start));
       if info.iter > maxiter
           break;
       end
   end
   t_end = clock();
   info.time = etime(t_end, t_start);
   if info.iter < maxiter
       info.solved = 1;
   else
       info.solved = 0;
   end
end