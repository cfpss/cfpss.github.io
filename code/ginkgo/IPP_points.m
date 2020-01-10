function S = IPP_points()

% Example 2.  IPP
% Generate data from 2-d space, let lambda(x,y)=exp(-(x-9).^2./10-(y-9).^2./10);
clear;
addpath('C:\Users\hongxiao\Documents\teaching\stat6474_2017\lectures\Lecture_22_23_Gaussian_process_regression\Lec2_GP_application\new_organized_code')
s = RandStream('mt19937ar','Seed',152);
RandStream.setGlobalStream(s);

ngrid=[30,20]; % num of grid points on each dimension.
maxTx=20;
maxTy=20;
lam_max=0.01; % put 0.1 for 5 trees maximum of the lam(s)
[t1,t2] = meshgrid(linspace(0,maxTx,ngrid(1)), linspace(0,maxTy,ngrid(2)));     

figure()
lambda_true = lam_max*exp(-(t1-9).^2/22-(t2-5).^2/22)+lam_max*exp(-(t1-22).^2/15-(t2-15).^2/22);

figure()
surf(t1,t2,lambda_true)
title('True lambda(s) surface');
xlabel('x')
ylabel('y')

sig_fn = inline('exp(-(x(:,1)-9).^2/22-(x(:,2)-5).^2/22) + exp(-(x(:,1)-15).^2/22-(x(:,2)-15).^2/22)'); 
T=[0,maxTx;0,maxTy];
[A,vol_T]=inhomo_poisson_process_sampler(lam_max,T,sig_fn);


imagesc(linspace(0,maxTx,ngrid(1))',linspace(0,maxTy,ngrid(2))',lambda_true)
hold on
plot(A(:,1),A(:,2),'.k','MarkerSize',20)    
hold off
title('IPP');

S=A;

end