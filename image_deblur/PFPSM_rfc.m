function out = PFPSM_rfc(H,F,mu,opts)
% out = PFPSM_rf(H,F,mu,opts)
%
% Alternating Directions Method (ADM) applied to TV/L2.
%
% Suppose the data accuquisition model is given by: F = K*Xbar + Noise,
% where Xbar is an original image, K is a convolution matrix, Noise is
% additive noise, and F is a blurry and noisy observation. To recover
% Xbar from F and K, we solve TV/L2 (ROF) model
%
% ***     min_X \sum_i ||Di*X|| + mu/2*||K*X - F||^2      ***
%
% Inputs:
%         H  ---  convolution kernel representing K
%         F  ---  blurry and noisy observation
%         mu ---  model prameter (must be provided by user)
%         opts --- a structure containing algorithm parameters {default}
%                 * opst.beta    : a positive constant {10}
%                 * opst.gamma   : a constant in (0,1.618] {1.618}
%                 * opst.maxitr  : maximum iteration number {500}
%                 * opst.relchg  : a small positive parameter which controls
%                                  stopping rule of the code. When the
%                                  relative change of X is less than
%                                  opts.relchg, then the code stops. {1.e-3}
%                 * opts.print   : print inter results or not {1}
%
% Outputs:
%         out --- a structure contains the following fields
%                * out.snr   : SNR values at each iteration
%                * out.f     : function valuse at each itertion
%                * out.relchg: the history of relative change in X
%                * out.sol   : numerical solution obtained by this code
%                * out.itr   : number of iterations used
%

% Copyright (c), May, 2009
%       Junfeng Yang, Dept. Math., Nanjing Univiversity
%       Wotao Yin,    Dept. CAAM, Rice University
%       Yin Zhang,    Dept. CAAM, Rice University


% PFPSM_rf adopts a slightly conservative strategy for the parameter $\alpha_{k}$, 
% i.e., we set $\alpha_{k}= \frac{1}{2m\tau+m+1}$ for all $k\geq 0$.

[m,n,d3] = size(F);

if d3 == 3
    error('Error, Grayscale image only!');
end

if nargin < 4; opts = []; end
opts = getopts(opts);

C = getC;
[D,Dt] = defDDt;

% initialization
X = F;
Lam1 = zeros(m,n);
Lam2 = Lam1;
beta = opts.beta;
gamma = opts.gamma;
tau = opts.tau;
print = opts.print;

% finite diff
[D1X,D2X] = D(X);
Y1 = D1X;
Y2 = D2X;
f = fval;

out.snr = [];
out.relchg = [];
out.f = f;

%% Main loop
for ii = 1:opts.maxitr
    %% Predictor
    % ==================
    %   Shrinkage Step
    % ==================
%     Sum1p = Y1-D1X;
%     Sum2p = Y2-D2X;
    Z1 = D1X + Lam1/beta;
    Z2 = D2X + Lam2/beta;
    V = Z1.^2 + Z2.^2;
    V = sqrt(V);
    V(V==0) = 1;
    V = max(V - 1/beta, 0)./V;
    barY1 = Z1.*V; 
    barY2 = Z2.*V; 
    
    % ==================
    %     X-subprolem
    % ==================
    Xp = X;
    X = (mu*C.KtF - Dt(Lam1,Lam2))/beta + Dt(Y1,Y2); 
    X = fft2(X)./(C.eigsDtD + (mu/beta)*C.eigsKtK);
    barX = real(ifft2(X));
    
    % finite diff.
    [barD1X,barD2X] = D(barX);
    Sum1 = barY1 - barD1X;
    Sum2 = barY2 - barD2X;
    
    % ==================
    %    Update Lam
    % ==================
    barLam1 = Lam1 - tau*beta*Sum1;
    barLam2 = Lam2 - tau*beta*Sum2;
    
    % Compute correction steplength
%     tpd = norm(barY1-Y1,'fro')^2+norm(barY2-Y2,'fro')^2+ ... 
%         norm(barD1X-D1X,'fro')^2+norm(barD2X-D2X,'fro')^2; % the distance of two point difference
%     varphi = tpd+norm(Sum1p,'fro')^2+norm(Sum2p,'fro')^2;
%     phi = tpd+tau*(norm(Sum1,'fro')^2+norm(Sum2,'fro')^2)+ ... 
%         norm(Sum1-Sum1p,'fro')^2+norm(Sum2-Sum2p,'fro')^2;
%     alphak = gamma*varphi/phi;
    alphak = 1/(2*2*tau+3);  % alphak = 1/(2*m*tau+m+1);
    
    % Correction step
    X = (1-alphak)*Xp+alphak*barX;
    
    snrX = snr(X);
    out.snr = [out.snr; snrX];
    relchg = norm(X - Xp,'fro')^2/norm(Xp,'fro')^2;
    out.relchg = [out.relchg; relchg];
    
    if print
        fprintf('Iter: %d, snrX: %4.2f, relchg: %4.2e\n',ii,snrX,relchg);
    end
    
    % ====================
    % Check stopping rule
    % ====================
    if relchg < opts.relchg
        out.sol = X;
        out.itr = ii;
        [D1X,D2X] = D(X);
        f = fval;
        out.f = [out.f; f];
        return
    end
    
    % finite diff.
    [D1X,D2X] = D(X);
    
    f = fval;
    out.f = [out.f; f];
    
    % ==========================
    %    Update Y1,Y2,Lam1,Lam2
    % ==========================
    Y1 = (1-alphak)*Y1+alphak*barY1;
    Y2 = (1-alphak)*Y2+alphak*barY2;
    Lam1 = (1-alphak)*Lam1+alphak*barLam1;
    Lam2 = (1-alphak)*Lam2+alphak*barLam2;
    
end
out.sol = X;
out.itr = ii;
out.exit = 'Exist Normally';
if ii == opts.maxitr
    out.exit = 'Maximum iteration reached!';
end

    function opts = getopts(opts)
        
        if ~isfield(opts,'maxitr')
            opts.maxitr = 500;
        end
        if ~isfield(opts,'beta')
            opts.beta = 10;
        end
        if ~isfield(opts,'gamma')
            opts.gamma = 1;  
        end
        if ~isfield(opts,'tau')
            opts.tau = 0.1;  % 0.97 for Circles and Barbara; 0.7 for Cameraman; 0.1 for Chart;
        end
        if ~isfield(opts,'relchg')
            opts.relchg = 1.e-5;
        end
        if ~isfield(opts,'print')
            opts.print = 0;
        end
    end


%% nested functions

    function C = getC
        sizeF = size(F);
        C.eigsK = psf2otf(H,sizeF);
        C.KtF = real(ifft2(conj(C.eigsK) .* fft2(F))); 
        C.eigsDtD = abs(psf2otf([1,-1],sizeF)).^2 + abs(psf2otf([1;-1],sizeF)).^2;
        C.eigsKtK = abs(C.eigsK).^2;
    end

    function f = fval
        f = sum(sum(sqrt(D1X.^2 + D2X.^2))); %sum all elements of sqrt(D1X.^2 + D2X.^2)
        KXF = real(ifft2(C.eigsK .* fft2(X))) - F;
        f = f + mu/2 * norm(KXF,'fro')^2;
    end

    function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);
        
        function [Dux,Duy] = ForwardD(U)
            % Forward finite difference operator
            Dux = [diff(U,1,2), U(:,1) - U(:,end)];
            Duy = [diff(U,1,1); U(1,:) - U(end,:)];
        end
        
        function DtXY = Dive(X,Y)
            % Transpose of the forward finite difference operator
            DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
        end
    end

end