% Solve the exchange probelm:
%
%   Minimize    f_1(x_1)+ ... + f_N(x_N)
%   subject to  x_1 + ... + X_N = 0
%
% where f_i(x_i)=0.5*||C_i*x_i-d_i||^2+xi*||x_i||_1.
%-----------------------------------------------------
% This MATLAB script is adapted from the following paper:
% Deng, W., Lai, M.J., Peng, Z., Yin, W.: Parallel multi-block ADMM with 
% O(1/k) convergence. J. Sci. Comput. 71(2), 712-736 (2017)
% see the authors' homepage or https://github.com/ZhiminPeng/Jacobi-ADMM

clear;clc

%seed = 2014; % use fixed seed
seed = sum(100*clock); % use clock seed
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

% Problem size
n = 100;        % length of x_i
m = 80;         % length of d_i
N = 100;        % number of x_i's
tol = -1;       % tolerance
maxit = 1000;   % max number of iterations
nrun = 10;      % number of runs
xi = 5e-4;      % the model parameter 

% Record residuals
Res_PFPSM_rf = zeros(maxit,nrun);
Res_PPSM = zeros(maxit,nrun);

% Record objective values
Obj_PFPSM_rf = zeros(maxit,nrun);
Obj_PPSM = zeros(maxit,nrun);

% Run test
for j = 1:nrun 
    fprintf('----- test %4i -----\n', j);
    %% Generate data (C,x,d)
    X0 = full(sprandn(n,N,0.5));     % X0 = randn(n,N);
    X0(:,N) = -sum(X0(:,1:N-1),2);
    C = cell(N,1);
    d = cell(N,1);
    CtC = cell(N,1);
    Ctd = cell(N,1);
    s = cell(N,1);       % s_i*I_n-rho C_i'*C_i, s_i¡Ý||C_i||^2
    for i = 1:N
        C{i} = randn(m,n);      
        d{i} = C{i}*X0(:,i);
        CtC{i} = C{i}'*C{i};
        Ctd{i} = C{i}'*d{i};
        s{i} = norm(C{i})^2;
    end
    %% PFPSM_rf
    opts1.rho = 3;
    opts1.tau = 1.1;
    opts1.gamma = 1;
    opts1.maxit = maxit;
    opts1.tol = tol;
    opts1.xi = xi;
    [X,~,Out1] = ExchgQuad_PFPSM_rf(C,d,CtC,Ctd,s,opts1);
    err = norm(X-X0,'fro')/norm(X0,'fro');
    fprintf('PFPSM_rf: iter = %4i, relative error = %e\n',...
        Out1.iter,err)
    Obj_PFPSM_rf(:,j) = Out1.objValue;
    Res_PFPSM_rf(:,j) = Out1.residual;
    
    %% PPSM
    opts2.rho = 1;
    opts2.gamma = 1;
    opts2.maxit = maxit;
    opts2.tol = tol;
    opts2.xi = xi;
    [X,~,Out2] = ExchgQuad_PPSM(C,d,CtC,Ctd,s,opts2);
    err = norm(X-X0,'fro')/norm(X0,'fro');
    fprintf('PPSM    : iter = %4i, relative error = %e\n',...
        Out2.iter,err)
    Obj_PPSM(:,j) = Out2.objValue;
    Res_PPSM(:,j) = Out2.residual;
    
%     %% Plot results
%     figure(1);
%     % Plot objective values
%     subplot(1,2,1), semilogy(1:Out1.iter, Out1.objValue,'b-');hold on
%     semilogy(1:Out2.iter, Out2.objValue,'k-.');
% %     semilogy(1:Out3.iter, Out3.objValue,'m-.'); 
% %     semilogy(1:Out4.iter, Out4.objValue,'g-');
%     hold off
%     xlabel('Iteration');
%     ylabel('Objective Value');
%     legend('PFPSM\_rf','PPSM')
%     % Plot residuals
%     subplot(1,2,2), semilogy(1:Out1.iter, Out1.residual,'b-');
%     hold on    
%     semilogy(1:Out2.iter, Out2.residual,'k-.');
% %     semilogy(1:Out3.iter, Out3.residual,'m-.');
%     %semilogy(1:Out4.iter, Out4.residual,'g-');
%     hold off
%     xlabel('Iteration');
%     ylabel('Residual');
%     legend('PFPSM\_rf','PPSM')
end
%return
%% Plot average results
figure(1);
lw = 2; % set line width
% Plot objective values
% subplot(1,2,1);
semilogy(1:maxit, geomean(Obj_PFPSM_rf,2),'b-','LineWidth',lw);
hold on
semilogy(1:maxit, geomean(Obj_PPSM,2),'k-.','LineWidth',lw);
% semilogy(1:maxit, geomean(Obj_CorrJADMM,2),'m-.','LineWidth',lw);
hold off
xlabel('Iteration','FontSize',12);
ylabel('Objective Value','FontSize',12);
legend('PFPSM\_rf','PPSM')

% Plot residuals
figure(2);
% subplot(1,2,2);
semilogy(1:maxit, geomean(Res_PFPSM_rf,2),'b-','LineWidth',lw);
hold on
semilogy(1:maxit, geomean(Res_PPSM,2),'k-.','LineWidth',lw);
% semilogy(1:maxit, geomean(Res_CorrJADMM,2),'m-.','LineWidth',lw); 
hold off
xlabel('Iteration','FontSize',12);
ylabel('Residual','FontSize',12);
legend('PFPSM\_rf','PPSM')

% Save data
clear X0 C d X Out X2 Out2 err;
%save ExchgQuad.mat