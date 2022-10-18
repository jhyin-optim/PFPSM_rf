% use PFPSM_rf proposed by Yin et al. to solve the following problem
% min f(x):=||x||_1, s.t. Ax=b, lb¡Üx¡Üub,
% where A¡ÊR^{m¡Án}. 
% For this 

function out = PFPSMrf(xs,lb,ub,tau,para)

global A b m n

% start clock
tic

% Set parameters
% tau = para.tau;
gamma = para.gamma;
beta = para.beta;
Itrmax = para.Itrmax;
epsilon = para.epsilon;

% Set initial point
x0 = zeros(n,1);
lam0 = zeros(m,1);

% define function handle
Proj = @(x,y,z) max(x,min(y,z));
pri_feas = @(x) A*x-b;

% save some useful product
diag_AtA = diag(A'*A);       % diag(X) is the main diagonal of X. 
Diag_AtA = diag(diag_AtA);   % diag(diag(X)) is a diagonal matrix.

for k=1:Itrmax
    
    % solve x-subproblems
%     hat_x1 = zeros(size(x0));
%     for i=1:n
%         p = x0(i)-A(:,i)'*(pri_feas(x0)-lam0/beta)/(AtA(i));
%         mu = 1/(beta*AtA(i));
%         hat_x1(i) = Proj(lb(i),ub(i),p-Proj(-mu,mu,p));
%     end
    P = x0-Diag_AtA\A'*(pri_feas(x0)-lam0/beta);
    hat_x1 = Proj(lb,ub,P-Proj(-1./(beta*diag_AtA),1./(beta*diag_AtA),P));

    % check the stopping condition
    if max(norm(hat_x1-x0,inf),norm(pri_feas(hat_x1),inf))<epsilon
        break;
    end
    
    % update multiplier
    hat_lam1 = lam0-tau*beta*pri_feas(hat_x1);
    
    % generate the next iterate
    sum = beta*norm(A*diag(hat_x1-x0),'fro')^2;
%     sum = 0;
%     for i=1:n
%         sum = sum+beta*norm(A(:,i)*(hat_x1(i)-x0(i)))^2;
%     end
    varphi = sum+beta*norm(pri_feas(x0))^2;
    psi = tau*beta*norm(pri_feas(hat_x1))^2+sum+beta*norm(A*(x0-hat_x1))^2;
    x0 = x0-gamma*varphi/psi*(x0-hat_x1);
    lam0 = lam0-gamma*varphi/psi*(lam0-hat_lam1);
    

end
out.Itr = k;
out.Tcpu = toc;
out.consm = norm(pri_feas(x0),inf);
out.obj = abs(norm(x0,1)-norm(xs,1));

