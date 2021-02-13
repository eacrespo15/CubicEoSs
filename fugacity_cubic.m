function [fug,rho] = fugacity_cubic(NC,T,P,zfeed,phasetype,solvertype,INITDATA)
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%This routine is just to make more explicit how to call the cubic EoSs to
%calculate fugacity coefficients.
[~,Z,fug,~,~,~,~] = cubic_eos(NC,T,P,zfeed,phasetype,solvertype,INITDATA.SQTC,INITDATA.kij,INITDATA.delta1,INITDATA.delta2,INITDATA.bci,INITDATA.aci,INITDATA.mfunc);
PINV=1/P;
V=Z*T*PINV; %This is V/R 
rho=(1/V)/1000; %Density in mol/L?
end

