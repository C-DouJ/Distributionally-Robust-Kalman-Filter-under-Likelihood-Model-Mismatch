%%%%Set up
add_path
clear all;
close all;
clc;
format long;

%%%%random state
randomselect='random'; % 'random' or 'repeatlast';
switch randomselect
  case 'repeatlast'
    load randstates;
    rand('state',randstate);
    randn('state',randnstate);
  case 'random'
    % do nothing
  otherwise
    rand('state',randomselect);
    randn('state',randomselect);
end
% %
randstate=rand('state');
randnstate=randn('state');
save randstates.mat randstate randnstate;

%%%%Model parameter
nxp=200; %Monte Carlo counts
nx=4;
nz=2;
T=1;        
q=1;
r=100;
F=[eye(2) T*eye(2);zeros(2) eye(2)];
H=[eye(2) zeros(2)];
% H=eye(4);
Q0=[T^3/3*eye(2) T^2/2*eye(2);T^2/2*eye(2) T*eye(2)]*q;
R0=r*eye(nz);
ts=200;

%%%%The choices of outlier parameters
pOutlier = 0.1;
U = 100;


t_KF = 0;
t_DRKF = 0;

for expt = 1:nxp
    
    fprintf('MC Run in Process = %d\n',expt); 
    
    %%%%Initial values
    x=[0;0;10;10];                          
    P=diag([10000 10000 100 100]);        
    Skk=utchol(P);                        

    %%%%KF (Kalman filter)
    xf=x+Skk*randn(nx,1);                 
    Pf=P;
    
    %%%%The Proposed Robust filter
    xDRKF = xf;
    PDRKF = Pf;

    %%%%Save data
    xA=x;
    xfA=xf;
    xDRKFA=xDRKF;

    PfA = Pf;
    PDRKFA = PDRKF;
    %%%%
    
    for t=1:ts
        %%%%Simulate true state and measurement
        SQ=utchol(Q0);  
        x=F*x+SQ*randn(nx,1);

        %%%%outlier
        if rand>pOutlier
            R=R0;
        else
            R=U*R0;
        end
        SR=utchol(R); 
        %random rotation perturbation
         H = [1 0 0 0;0 1 0 0];
         perturbation = 2;%deg
         if t >80 && t<120
            theta1 = perturbation/57.3*randn;
            H = [cos(theta1) -sin(theta1);sin(theta1) cos(theta1)]*[1 0 0 0;0 1 0 0] ;
         end

         z=H*x+SR*(randn(nz,1));

         %uniform distribution measurement nosie
%           z=H*x+ 400*rand(nz,1)-200; 
        %%%%Filtering
         H = [1 0 0 0;0 1 0 0];

        tic
        [xf,Pf,xf1,Pf1]=kf(xf,Pf,F,H,z,Q0,R0);
        t_KF = t_KF + toc/nxp/ts;

        tic
        [xDRKF,PDRKF] = DRKF(xDRKF,PDRKF,F,H,z,Q0,R0);
        t_DRKF = t_DRKF + toc/nxp/ts;

    
        %%%%Save data
        xA=[xA x];
        xfA=[xfA xf];
        xDRKFA=[xDRKFA xDRKF];


        PfA(:,:,t+1) = Pf;
        PDRKFA(:,:,t+1) = PDRKF;
%%     CalMSE
        mse_kf_1(1,t,expt)=(xA(1,t)-xfA(1,t))^2+(xA(2,t)-xfA(2,t))^2;
        mse_kf_2(1,t,expt)=(xA(3,t)-xfA(3,t))^2+(xA(4,t)-xfA(4,t))^2;


        mse_DRKF_1(1,t,expt)=(xA(1,t)-xDRKFA(1,t))^2+(xA(2,t)-xDRKFA(2,t))^2;
        mse_DRKF_2(1,t,expt)=(xA(3,t)-xDRKFA(3,t))^2+(xA(4,t)-xDRKFA(4,t))^2;

        
    end
end

%%%%%%%%%RMSE calculation
rmse_kf_1=sqrt(mean(mse_kf_1,3));
rmse_kf_2=sqrt(mean(mse_kf_2,3));

rmse_DRKF_1=sqrt(mean(mse_DRKF_1,3));
rmse_DRKF_2=sqrt(mean(mse_DRKF_2,3));

t_start = 80;
t_end = 140;
Armse_kf_1=(mean(rmse_kf_1(t_start:t_end)))
Armse_kf_2=(mean(rmse_kf_2(t_start:t_end)))

Armse_DRKF_1=(mean(rmse_DRKF_1(t_start:t_end)))
Armse_DRKF_2=(mean(rmse_DRKF_2(t_start:t_end)))


disp(["t_KF = "+num2str(t_KF)])
disp(["t_DRKF = "+num2str(t_DRKF)])

%%%%%%%RMSE curve
plot_figure
