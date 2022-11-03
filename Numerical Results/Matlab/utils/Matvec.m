%%********************************************************************
%% normytmp = norm_yinput(idx); 
%% Dsub = yinput(:,idx)*spdiags(1./normytmp',0,len,len);
%% alp =  weight(idx)./(par.sigma*normytmp); 
%% Copyright (c) 2021 by
%% Defeng Sun, Kim-Chuan Toh and Yancheng Yuan 
%% Citation:
%% [1] D.F. Sun, K.C. Toh, and Y.C. Yuan, Convex clustering: model, theoretical guarantee and efficient algorithm, Journal of Machine Learning Research, 22(9):1?32, 2021.
%% [2] Y.C. Yuan, D.F. Sun, and K.C. Toh, An efficient semismooth Newton based algorithm for convex clustering, International Conference on Machine Learning (ICML) 2018.
%%********************************************************************
function My = Matvec(y,par,Ainput)

    idx = par.nzidx;     
    len = length(idx);
    Ay = Ainput.Amap(y);    
    if (len > 0)
       Aytmp = Ay(:,idx);        
       alpha  = par.alpha;
       Dsub = par.Dsub; 
       rho = alpha.*sum(Aytmp.*Dsub);      
       Ay(:,idx) = bsxfun(@times,Aytmp,alpha)-bsxfun(@times,Dsub,rho);
    end
    My = y + par.sigma*Ainput.ATmap(Ay);    
%%********************************************************************
