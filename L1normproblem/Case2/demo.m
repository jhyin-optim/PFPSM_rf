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

noise = 0; % If noise=1, then there is a Gauss white noise; otherwise, there is no noise

% Set parameters for PFPSM_rf
para1.tau = 1.1;
para1.gamma = 1;
para1.Itrmax = 1000;
para1.epsilon = 1e-6;

% Set parameters for FJDADM
para2.gamma = 1;
para2.Itrmax = 1000;
para2.epsilon = 1e-6;

% % Set parameters for ALBPS as the special case of PFPSM_rf with
% % C_i=0, tau=1
% para3.tau = 1;
% para3.gamma = 1;
% para3.Itrmax = 1000;
% para3.epsilon = 1e-6;

data = [10 100 1;20 200 2;50 500 5;100 1000 10;200 2000 20;500 5000 50;700 7000 70;900 9000 90];
[row,~] = size(data);

fid=fopen('mytext.txt','w');
for index=1:row
    m = data(index,1);
    n = data(index,2);
    s = data(index,3);  % s denotes the number of nonzero entries of x
    
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
        out1 = PFPSMrf(xs,para1);

        disp('Starting FJDADM')
        para2.beta = 10/sqrt(n);     % 10/sqrt(n)
        out2 = FJDADM(xs,para2);
        
%         disp('Starting ALBPS')
%         para3.beta = 10/sqrt(n);
%         out3 = ALBPS(xs,para3);
        
        progress_r=[progress_r; out1.Itr out1.Tcpu out1.obj out1.consm ... 
             out2.Itr out2.Tcpu out2.obj out2.consm];
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

