%*************************************************************************
% SSNAL: A semismooth Newton-CG augmented Lagrangian method for convex clustering
%
% Minor changes by Daniel Touw in order to prevent the code from printing
% intermediate results to reduce computation time.
%
% Copyright (c) 2021 by
% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
% Citation:
% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, theoretical guarantee and efficient algorithm, Journal of Machine Learning Research, 22(9):1?32, 2021.
% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton based algorithm for convex clustering, International Conference on Machine Learning (ICML) 2018.
%%*************************************************************************
function [obj,y,xi,z,info,runhist] = SNAL(Ainput,b,dim,weightVec,options,y,xi,z)

   rng('default');
   sigma = 1; 
   maxiter = 100;
   stoptol = 1e-6;
   printyes = 1;
   printminoryes=1;
   admm.iter = 0;
   admm.time = 0;
   admm.timecpu = 0;
   scale = 0;
   use_kkt = 1;
   if isfield(options,'maxiter');  maxiter  = options.maxiter; end
   if isfield(options,'stoptol');  stoptol  = options.stoptol; end
   if isfield(options,'printyes'); printyes = options.printyes; end
   if isfield(options,'printminoryes'); printminoryes = options.printminoryes; end
   if isfield(options,'sigma'); sigma = options.sigma; end
   if isfield(options,'use_kkt'); use_kkt = options.use_kkt; end
%%   
%% Amap and ATmap
%%
   Fnorm = @(x) mexFnorm(x);
   tstart = clock;
   tstart_cpu = cputime;
   if ~exist('z','var') || ~exist('xi','var') || ~exist('y','var')
      xi = b; z = sparse(dim.d,dim.E); y = Ainput.Amap(xi);
   end
   %% phase I
   if (options.admm_iter > 0)
       admmopt.stoptol = stoptol;
       admmopt.maxiter = options.admm_iter;
       admmopt.use_kkt = use_kkt;
       admmopt.solver = 'direct';
       
       % Prevent ADMM from printing information
       admmopt.printyes = 0;
       admmopt.printminoryes = 0;
       
       [obj,y,xi,z,info_admm,runhist_admm] = ADMM(Ainput,b,dim,weightVec,admmopt,y,xi,z);
       admm.iter = admm.iter + info_admm.iter;
       admm.time = admm.time + info_admm.time;
       admm.timecpu = admm.timecpu + info_admm.time_cpu;
       sigma = min(info_admm.sigma,10);
       if (info_admm.eta < stoptol)
          fprintf('\n Problem solved in Phase I \n');
          info = info_admm;
          info.dim = dim;
          info.minx = min(min(xi));
          info.maxx = max(max(xi));
          info.relgap = info_admm.relgap;
          info.iter = 0;
          info.time = admm.time;
          info.time_cpu = admm.timecpu;
          info.admmtime = admm.time;
          info.admmtime_cpu = admm.timecpu;
          info.admmiter = admm.iter;
          info.eta = info_admm.eta;
          info.etaorg = info_admm.etaorg;
          info.obj = obj;
          info.maxfeas = max([info_admm.dualfeasorg, info_admm.primfeasorg]);
          runhist = runhist_admm;
          return
       end
   end
   %% Initialization
   breakyes = 0;
   normy = Fnorm(y);
   Axi = Ainput.Amap(xi); Atz = Ainput.ATmap(z);
   Rp = Axi - y;
   proj_z = proj_l2(z,weightVec);
   Rd = z - proj_z;
   primfeas = Fnorm(Rp)/(1+normy);
   dualfeas = Fnorm(Rd)/(1 + Fnorm(z));
   maxfeas = max(primfeas, dualfeas);
   primfeasorg = primfeas;
   dualfeasorg = dualfeas;
   primobj = 0.5*Fnorm(xi-b)^2 + sum(weightVec.*(sqrt(sum(Axi.*Axi))));
   dualobj = -0.5*Fnorm(Atz)^2 + sum(sum(b.*Atz));
   relgap = abs(primobj - dualobj)/( 1+abs(primobj)+abs(dualobj));
   if printyes
        fprintf('\n *******************************************************');
        fprintf('******************************************');
        fprintf('\n \t\t   Phase II: NAL_Clustering ');
        fprintf('\n ******************************************************');
        fprintf('*******************************************\n');
        if printminoryes
           fprintf(' d = %3.0f, n = %3.0f',dim.d, dim.n);
           fprintf('\n ---------------------------------------------------');
        end
        fprintf('\n  iter|  [pinfeas  dinfeas pinforg  dinforg relgaporg]    primobj         daulobj   |');
        fprintf(' time | sigma |rankS|');
        fprintf('\n*****************************************************');
        fprintf('**************************************************************');
        fprintf('\n #%3.1d|  %3.2e %3.2e %3.2e %3.2e %- 3.2e %- 8.7e %- 8.7e  %5.1f',...
           0,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,primobj,dualobj,admm.time); 
        fprintf('  %3.2e ',sigma);
   end
   if (maxfeas < max(1e-6, stoptol))
        if use_kkt == 1
            grad = Atz + xi - b;
            eta = Fnorm(grad)/(1 + Fnorm(xi));
            ypz = y + z;
            [ypz_prox] = prox_l2(ypz, weightVec);
            res_vec = y - ypz_prox;
            eta = eta + Fnorm(res_vec)/(1+Fnorm(y));
            etaorg = eta;
        else
            primobj = 0.5*Fnorm(xi-b)^2 + sum(weightVec.*(sqrt(sum(Axi.*Axi))));
            dualobj = -0.5*Fnorm(Atz)^2 + sum(sum(b.*Atz));
            relgap = abs(primobj - dualobj)/( 1+abs(primobj)+abs(dualobj));
            eta = relgap;
            etaorg = eta;
        end
        if (eta < stoptol)
             breakyes = 1;
             msg = 'converged';
        end
   end
   if breakyes == 1
       obj = [primobj, dualobj];
       runhist.d = dim.d;
       runhist.n = dim.n;
       runhist.E = dim.E;
       ttime = etime(clock,tstart);
       ttime_cpu = cputime - tstart_cpu;
       runhist.iter = 0;
       runhist.totaltime = ttime;
       runhist.primobjorg = primobj; 
       runhist.dualobjorg = dualobj;
       runhist.maxfeas = max([dualfeasorg, primfeasorg]);
       runhist.eta = eta;
       runhist.etaorg = etaorg;
       info.d = dim.d;
       info.n = dim.n;
       info.E = dim.E;
       info.minx = min(min(z));
       info.maxx = max(max(z));
       info.relgap = relgap;
       info.iter = 0;
       info.time = ttime;
       info.time_cpu = ttime_cpu;
       info.admmtime = admm.time;
       info.admmtime_cpu = admm.timecpu;
       info.admmiter = admm.iter;
       info.eta = eta;
       info.etaorg = etaorg;
       info.obj = obj;
       info.sigma = sigma; 
       info.maxfeas = max([dualfeasorg, primfeasorg]);
       if (printminoryes) 
           if ~isempty(msg); fprintf('\n %s',msg); end
               fprintf('\n--------------------------------------------------------------');
               fprintf('------------------');
               fprintf('\n  admm iter = %3.0d, admm time = %3.1f', admm.iter, admm.time);
               fprintf('\n  number iter = %2.0d',0);       
               fprintf('\n  time = %3.2f',ttime);       
               fprintf('\n  time per iter = %5.4f',ttime); 
               fprintf('\n  cputime = %3.2f', ttime_cpu);
               fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj, relgap);       
               fprintf('\n  primfeasorg = %3.2e, dualfeasorg = %3.2e',...
               primfeasorg, dualfeasorg); 
               fprintf('\n  eta = %3.2e, etaorg = %3.2e', eta, etaorg);
               fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
               info.minx,info.maxx); 
               fprintf('\n--------------------------------------------------------------');
               fprintf('------------------\n');
       end
   end
   %% ssncg
   use_SSNCG = 1;
   if use_SSNCG
      parNCG.matvecfname = 'matvecIpxi';
      parNCG.sigma = sigma;
      parNCG.tolconst = 0.5;
      parNCG.dim = dim;
   end
   maxitersub = 10;
   prim_win = 0;
   dual_win = 0;
   ssncgopt.tol = stoptol;
   ssncgopt.precond = 0;
   ssncgopt.bscale = 1;
   ssncgopt.cscale = 1;
   for iter = 1:maxiter
      zold = z;      
      parNCG.sigma = sigma;
      if primfeas < 1e-5
         maxitersub = max(maxitersub,30);
      elseif primfeas < 1e-3
         maxitersub = max(maxitersub,30);
      elseif primfeas < 1e-1
         maxitersub = max(maxitersub,20);
      end
      ssncgopt.maxitersub = maxitersub;        
      [result, y,Axi,xi,parNCG,info_NCG] = evalc("SSNCG(b,Ainput,z,Axi,xi,weightVec,parNCG,ssncgopt)");
      %[y,Axi,xi,parNCG,info_NCG] = SSNCG(b,Ainput,z,Axi,xi,weightVec,parNCG,ssncgopt);
      if (info_NCG.breakyes < 0)
         parNCG.tolconst = max(parNCG.tolconst/1.06,1e-3);
      end
      Rp = Axi - y;
      z = zold + sigma*Rp;
      Atz = Ainput.ATmap(z);
      normRp = Fnorm(Rp); normy = Fnorm(y);
      primfeas = normRp/(1+normy);
      primfeasorg = primfeas;
      proj_z = proj_l2(z,weightVec);
      Rd = z - proj_z;
      dualfeas = Fnorm(Rd)/(1 + Fnorm(z));
      dualfeasorg = dualfeas;
      maxfeas = max(primfeas,dualfeas);
      maxfeasorg = maxfeas;
      runhist.dualfeas(iter+1) = dualfeas;
      runhist.primfeas(iter+1) = primfeas;
      runhist.maxfeas(iter+1)  = maxfeas;
      runhist.primfeasorg(iter) = primfeasorg;
      runhist.maxfeasorg(iter)  = maxfeasorg;
      runhist.sigma(iter) = sigma;
      runhist.rankS(iter) = sum(parNCG.rr);
      %%---------------------------------------------------------
%% check for termination
%%---------------------------------------------------------    
      if (maxfeas < max(1e-6, stoptol))
          if use_kkt == 1
            grad = Atz + xi - b;
            eta = Fnorm(grad)/(1 + Fnorm(xi));
            ypz = y + z;
            [ypz_prox] = prox_l2(ypz, weightVec);
            res_vec = y - ypz_prox;
            eta = eta + Fnorm(res_vec)/(1+Fnorm(y));
            etaorg = eta;
          else
            primobj = 0.5*Fnorm(xi-b)^2 + sum(weightVec.*(sqrt(sum(Axi.*Axi))));
            dualobj = -0.5*Fnorm(Atz)^2 + sum(sum(b.*Atz));
            relgap = abs(primobj - dualobj)/( 1+abs(primobj)+abs(dualobj));
            eta = relgap;
            etaorg = eta;
          end
          if (eta < stoptol)
             breakyes = 1;
             msg = 'converged';
          end
      end
          
%%--------------------------------------------------------    
%% print results
%%--------------------------------------------------------
      if printyes
         if maxfeas > stoptol
            primobj = 0.5*Fnorm(xi-b)^2 + sum(weightVec.*(sqrt(sum(Axi.*Axi))));
            dualobj = -0.5*Fnorm(Atz)^2 + sum(sum(b.*Atz));
            relgap = abs(primobj - dualobj)/( 1+abs(primobj)+abs(dualobj));
         end
         ttime = etime(clock,tstart);
         if (printyes)
            fprintf('\n %5.0d| [%3.2e %3.2e] [%3.2e %3.2e]  %- 3.2e| %- 5.4e %- 5.4e |',...
               iter,primfeas,dualfeas,primfeasorg,dualfeasorg,relgap,primobj,dualobj); 
            fprintf(' %5.1f| %3.2e|',ttime,sigma);                
            if (iter >= 1)
               fprintf('%3.0d|',sum(parNCG.rr));
            end
            if exist('eta'); fprintf('\n \t [ eta = %3.2e, etaorg = %3.2e]',eta,etaorg);end
         end
% 	     if (rem(iter,3*1)==1) 
%             normz = Fnorm(z); normAxi = Fnorm(Axi); normy = Fnorm(y);
%             if (printyes)
%                fprintf('\n        [normx,Atxi,y =%3.2e %3.2e %3.2e]',...
%                normx,normAxi,normy);
%             end
%          end
         runhist.primobj(iter)   = primobj;
         runhist.dualobj(iter)   = dualobj;
         runhist.time(iter)      = ttime; 
         runhist.relgap(iter)    = relgap;
      end    
     if (breakyes > 0) 
        fprintf('\n  breakyes = %3.1f, %s',breakyes,msg); 
        break; 
     end
     if (primfeasorg < dualfeasorg); 
         prim_win = prim_win+1; 
      else
         dual_win = dual_win+1; 
      end   
      if (iter < 10) 
         sigma_update_iter = 2;
      elseif iter < 20
         sigma_update_iter = 3;
      elseif iter < 200
         sigma_update_iter = 3;
      elseif iter < 500
         sigma_update_iter = 10;
      end
      sigmascale = 5;
      sigmamax = 1e5;
      if (rem(iter,sigma_update_iter)==0)
   	     sigmamin = 1e-4; 
	     if prim_win > max(1,1.2*dual_win)
            prim_win = 0;
            sigma = max(sigmamin,sigma/sigmascale);
         elseif dual_win > max(1,1.2*prim_win)
            dual_win = 0;
            sigma = min(sigmamax,sigma*sigmascale);
         end
      end
   end
%%-----------------------------------------------------------------
%% recover orignal variables
%%-----------------------------------------------------------------
   if (iter == maxiter)
      msg = ' maximum iteration reached';
      primobj = 0.5*Fnorm(xi-b)^2 + sum(weightVec.*(sqrt(sum(Axi.*Axi))));
      dualobj = -0.5*Fnorm(Atx)^2 + sum(sum(b.*Atx));
      relgap = (primobj-dualobj)/( 1+abs(primobj)+abs(dualobj));
      if use_kkt == 1
          grad = Atz + xi - b;
          eta = Fnorm(grad)/(1 + Fnorm(xi));
          ypz = y + z;
          [ypz_prox] = prox_l2(ypz, weightVec);
          res_vec = y - ypz_prox;
          eta = eta + Fnorm(res_vec)/(1+Fnorm(y));
          etaorg = eta;
      else
        eta = relgap;
        etaorg = eta;
      end
      info.termcode = 3;
   end
   obj = [primobj, dualobj];
   runhist.d = dim.d;
   runhist.n = dim.n;
   runhist.E = dim.E;
   ttime = etime(clock,tstart);
   ttime_cpu = cputime - tstart_cpu;
   runhist.iter = iter;
   runhist.totaltime = ttime;
   runhist.primobjorg = primobj; 
   runhist.dualobjorg = dualobj;
   runhist.maxfeas = max([dualfeasorg, primfeasorg]);
   runhist.eta = eta;
   runhist.etaorg = etaorg;
   info.d = dim.d;
   info.n = dim.n;
   info.E = dim.E;
   info.minx = min(min(z));
   info.maxx = max(max(z));
   info.relgap = relgap;
   info.iter = iter;
   info.time = ttime;
   info.time_cpu = ttime_cpu;
   info.admmtime = admm.time;
   info.admmtime_cpu = admm.timecpu;
   info.admmiter = admm.iter;
   info.eta = eta;
   info.etaorg = etaorg;
   info.obj = obj;
   info.sigma = sigma; 
   info.maxfeas = max([dualfeasorg, primfeasorg]);
   if (printminoryes) 
      if ~isempty(msg); fprintf('\n %s',msg); end
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------');
      fprintf('\n  admm iter = %3.0d, admm time = %3.1f', admm.iter, admm.time);
      fprintf('\n  number iter = %2.0d',iter);       
      fprintf('\n  time = %3.2f',ttime);       
      fprintf('\n  time per iter = %5.4f',ttime/iter); 
      fprintf('\n  cputime = %3.2f', ttime_cpu);
      fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj, relgap);       
      fprintf('\n  primfeasorg = %3.2e, dualfeasorg = %3.2e',...
	      primfeasorg, dualfeasorg); 
      fprintf('\n  eta = %3.2e, etaorg = %3.2e', eta, etaorg);
      fprintf('\n  min(X)    = %3.2e, max(X)    = %3.2e',...
          info.minx,info.maxx); 
      fprintf('\n--------------------------------------------------------------');
      fprintf('------------------\n');
   end
%%**********************************************************************

 
   