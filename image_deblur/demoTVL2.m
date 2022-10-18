%
% Demo TV/L2 solve
%
% Suppose the data accuquisition model is given by: F = K*Xbar + Noise,
% where Xbar is an original image, K is a convolution matrix, Noise is
% additive noise, and F is a blurry and noisy observation. To recover
% Xbar from F and K, we solve TV/L2 (ROF) model
%
% ***     min_X \sum_i ||Di*X|| + mu/2*||K*X - F||^2      ***
%

clear; close all;
% path(path,genpath(pwd));
addpath('./Images')
rng(2016);

%% generate data -- see nested function below
Ima = {'cameraman.tif','barbara.png','chart.tiff','circles.tif', ... 
            'housergb.png','lena.png','shape.jpg','TwoCircles.tif'};
index = 3;
I = double(imread(Ima{index}))/255;   
        
[m,n] = size(I);
% H = fspecial('disk',7);
H = fspecial('average',10);  % 此时Ima{6}效果不好，12对应效果好 
% H = fspecial('average',15); %创建一个15×15的均值滤波器 
        
sigma = 1.e-2;     %sigma = 1.5e-3 % sigma = 1.e-3;
Bn = imfilter(I,H,'circular','conv') + sigma*randn(m,n);

mu = 5.e3;   % mu = 5.e4;
snrBn = snr(Bn,I);  % The initial value of SNR 
%
fprintf(' ********** Comparison starts **********\n');

%% Run PFPSM_rf
fprintf(' ********** run PFPSM_rf ***************\n');
t1 = cputime;
out1 = PFPSM_rf(H,Bn,mu); 
t1 = cputime - t1;

%% Run FJADMM
fprintf(' ********** run FJADMM ***************\n');
t2 = cputime;
out2 = FJADMM(H,Bn,mu); 
t2 = cputime - t2;

%% Run PPSM
fprintf(' ********** run PPSM ***************\n');
t3 = cputime;
out3 = PPSM(H,Bn,mu); 
t3 = cputime - t3;

%% Run PPSM
fprintf(' ********** run PFPSM_rfc ***************\n');
t4 = cputime;
out4 = PFPSM_rfc(H,Bn,mu); 
t4 = cputime - t4;

% output results
fid=fopen('mytext.txt','w');
fprintf(fid,'%s & %.2f & %d/%.2f/%.2f & %d/%.2f/%.2f & %d/%.2f/%.2f & %d/%.2f/%.2f\\\\\n', ... 
        Ima{index},snrBn,out1.itr,t1,snr(out1.sol),out2.itr,t2,snr(out2.sol),out3.itr,t3,snr(out3.sol),out4.itr,t4,snr(out4.sol));   % mean(progress_r)
fclose(fid);
fprintf('PFPSM_rf SNR(Bn) %4.2fdB, SNR(Recovered) %4.2fdB,',snrBn,snr(out1.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t1,out1.itr);

fprintf('FJADMM SNR(Bn) %4.2fdB, SNR(Recovered) %4.2fdB,',snrBn,snr(out2.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t2,out2.itr);

fprintf('PPSM SNR(Bn) %4.2fdB, SNR(Recovered) %4.2fdB,',snrBn,snr(out3.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t3,out3.itr);

fprintf('PFPSM_rfc SNR(Bn) %4.2fdB, SNR(Recovered) %4.2fdB,',snrBn,snr(out4.sol));
fprintf(' CPU %4.2fs, Iteration %d\n\n',t4,out4.itr);

%% plot result
figure(1);
imshow(I,[]);
title('Original image');

figure(2);
imshow(Bn,[]);
title('Degraded image');

figure(3);
imshow(out1.sol,[]);
title('PFPSM\_rf');

figure(4);
imshow(out2.sol,[]);
title('FJADMM');

figure(5);
imshow(out2.sol,[]);
title('PPSM');

figure(6);
imshow(out4.sol,[]);
title('PFPSM\_rfc');

% figure(1);
% subplot(121); imshow(Bn,[]);
% title(sprintf('SNR %4.2fdB',snrBn),'fontsize',13);
% subplot(122); imshow(out1.sol,[]);
% title(sprintf('SNR %4.2fdB, CPU %4.2fs, It: %d',snr(out1.sol),t1,out1.itr),'fontsize',13);

% figure(2);
% subplot(121); plot(1:length(out1.snr),out1.snr,'b--',1:length(out1.snr),out1.snr,'r.');
% title('SNR history','fontsize',13);
% subplot(122); semilogy(1:length(out1.f)-1,out1.f(2:end),'b--',1:length(out1.f)-1,out1.f(2:end),'r.');
% title('Function values','fontsize',13);
% axis([0,inf,min(out1.f(2:end)),max(out1.f(2:end))]);
