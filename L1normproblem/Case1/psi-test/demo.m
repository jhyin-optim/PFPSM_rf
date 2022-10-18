% Consider the following l_1 norm model
% min f(x):=||x||_1, s.t. Ax=b, lb≤x≤ub,
% where A∈R^{m×n}.
% Use the FJDADM proposed by He, Hou and Yuan in 2015 
% and PFPSM_rf proposed by Yin et al. to solve the above problem.
%
% The test data see the following paper:
% He, B., Hou, L., Yuan, X.: On full Jacobian decomposition of the 
% augmented Lagrangian method for separable convex programming. 
% SIAM J. Optim. 25(4), 2274-2312 (2015)
%

clc
clear all
close all
% clf

% control the random number generator used by rand, randi, and randn.
rng(1,'twister');

global A b m n

noise = 0; % If noise=0, then there is no noise; otherwise, there is a Gauss white noise.

% Set parameters for PFPSM_rf
para1.tau = 1;
para1.gamma = 1;
para1.Itrmax = 1000;
para1.epsilon = 1e-6;

% % Set parameters for FJDADM
% para2.gamma = 1;
% para2.Itrmax = 1000;
% para2.epsilon = 1e-6;

% Set parameters for ALBPS as the special case of PFPSM_rf with
% C_i=0, tau=1
para3.tau = 1;
para3.gamma = 1;
para3.Itrmax = 1000;
para3.epsilon = 1e-6;

data = [10 25 2;20 50 5;50 100 10;100 300 15;200 500 20;500 1000 30; ... 
    1000 2000 50;2000 5000 100];
[row,~] = size(data);

fid=fopen('mytext.txt','w');
for index=1:row
    m = data(index,1);
    n = data(index,2);
    s = data(index,3);  % s denotes the number of nonzero entries of x
    % the lower bound
    lb = -ones(n,1);
    % the upper bound
    ub = ones(n,1);
    
    progress_r = [];
    for repeats=1:10
        % generate test data
        xs = zeros(n,1);  % initial original signal
        q = randperm(n);
        xs(q(1:s)) = sign(randn(s,1));  % generate original sparse signal and the componentwises are 1 or -1.
        A = normr(randn(m,n));  % each row of A is normalized as a vector with a length of 1;

        if noise==1
            b = A*xs+0.01*randn(m,1);
        else
            b = A*xs;
        end
        
%         % compute the Lipschitz constant
%         tic
%         Lipc = norm(A*A');
%         CPUtime = toc;
        
        disp('Starting PFPSM_rf')
        para1.beta = 3;
        out1 = PFPSMrf(xs,lb,ub,para1);

%         disp('Starting FJDADM')
%         para2.beta = 10/sqrt(n);     % 10/sqrt(n)
%         out2 = FJDADM(xs,lb,ub,para2);
        
        disp('Starting ALBPS')
        para3.beta = 3;
        out3 = ALBPS(xs,lb,ub,para3);
        
        progress_r=[progress_r; out1.Itr out1.Tcpu out1.obj out1.consm ... 
             out3.Itr out3.Tcpu out3.obj out3.consm];
    end
    fprintf(fid,'%d & %d & %d & %.0f/%.3f/%.3e/%.3e & %.0f/%.3f/%.3e/%.3e\\\\\n', ... 
        m,n,s,mean(progress_r));   % mean(progress_r)
end
fclose(fid);

% %==========================================================================
% %% plot
% figure(1)
% semilogy(0:para1.Itrmax,out1.objfunval,0:para2.Itrmax,out2.objfunval)
% hold on
% semilogy(0:para3.Itrmax,out3.objfunval)
% set(gca,'XTick',[0:200:1000]);%设置要显示坐标刻度
% title('目标函数变化趋势')
% legend('FISTA\_S','PG\_S','MFISTAS','Location','best') 
% xlabel('k')
% ylabel('f(x^k)-f\_min')
% hold off
% 
% figure(2)
% semilogy(0:para1.Itrmax,out1.distval,0:para2.Itrmax,out2.distval)
% hold on
% semilogy(0:para3.Itrmax,out3.distval)
% set(gca,'XTick',[0:200:1000]);%设置要显示坐标刻度
% title('距离变化趋势')
% legend('FISTA\_S','PG\_S','MFISTAS','Location','best') 
% xlabel('k')
% ylabel('||x^k-x\_s||')
% hold off

