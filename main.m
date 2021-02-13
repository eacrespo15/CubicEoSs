%Main file for calculations with SRK/PR cubic EoSs
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%This file provides two examples of how to use the SRK/PR codes


%% Example 1 -
Input required for calculations using Cubic EoSs
NC=3;                   %Number of components
indexes=[1 14 15];      %Index of the mixture components in the compounds database (Check get_module_database.m)
T=240;                  %Temperature (K)
P=5;                    %Pressure (MPa)
nfeed=[50 10 40];       %Feed (number of moles in the feed)
EoS=1;                  %EoS Model (0)SRK (1)Peng-Robinson
phasetype=1;            %Desired root: (1)Liquid (-1)Vapor (0)Minimum Gibbs energy
solvertype=1;           %Solver type: (1)Mollerup (2)Mollerups with root sort

%Convert the Feed into a normalized composition for n=1
zfeed=nfeed/sum(nfeed);
%Read the database and set up the correct cubic EoS
[SQTC,kij,delta1,delta2,bci,aci,mfunc]=init_cubic(NC,EoS,indexes);
%Carry the calculations
[phasefound,Z,fug,fugT,fugP,fugX,AUX]=cubic_eos(NC,T,P,zfeed,phasetype,solvertype,SQTC,kij,delta1,delta2,bci,aci,mfunc);

%% Particular notes of the developer
% NC=5;
% indexes=[1,2,3,14,15];
% T=204;
% P=4;
% z=[0.66 0.03 0.01 0.05 0.25];
% EoS=0;
% flashphase=0;   %Fugacities are always calculated for the phase with lower G for the trial composition
% solver_cubic=1; %Define solver type for cubic EoS
% 
% %Read the database
% [data]=get_module_database();
% %Get the critical properties
% Tc=data.TC(indexes);
% Pc=data.PC(indexes);
% W=data.ACENTRIC(indexes);

%Read the database and set up the correct cubic EoS
% [SQTC,kij,delta1,delta2,bci,aci,mfunc]=init_cubic(NC,EoS,indexes);
% INITDATA.SQTC=SQTC;
% INITDATA.kij=kij;
% INITDATA.delta1=delta1;
% INITDATA.delta2=delta2;
% INITDATA.bci=bci;
% INITDATA.aci=aci;
% INITDATA.mfunc=mfunc;
% 
% NF=3;                           %Assumed number of phases
% betas_guess=1/NF*ones(NF,1);    %Construct initial vector of beta
% %Define initial K-values
% K=zeros(1,NC);
% for i=1:NC
%    K(i)=Pc(i)/P*exp(5.373*(1+W(i))*(1-Tc(i)/T));
% end

%Initialize the fugacities
% LK_L1=log(K);
% LK_L1(end)=LK_L1(end)+1;
% LK_L2=log(K);
% LK_L2(1)=LK_L2(1)+1;
% FUG_V=ones(1,NC);
% FUG_L1=exp(LK_L1);
% FUG_L2=exp(LK_L2);
% FUGs=[FUG_V;FUG_L1;FUG_L2];

%Collect initial guess from rachford-rice solver
% [Betas,y]=MPh_Rachford_Rice(NF,z,FUGs,betas_guess);
% [Betas,y]=MF_Rachford_Rice(NF,z,FUGs,betas_guess);

%Apply successive substitutuion to solve the problem
%[Final_Betas,XPhase,rho_phase,niter,time]=MF_FLASH_MICHELSEN(P,T,z,y,Betas,flashphase,solver_cubic,INITDATA);


