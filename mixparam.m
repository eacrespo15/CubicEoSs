function [Bmix,Amix,Amix_dT,Amix_dT2,Amix_dX,Amix_dXdT] = mixparam(zfeed,kij,bii,ai0,ai1,ai2)
%This subroutines returns the A and B mixture parameters and the required
%composition and temperature derivatives
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Inputs:
%zfeed - Systems composition (1xNC)
%kij - Matrix of binary interaction parameters in the form of (1-kij)
%bii - bii values for every component (1xNC) 
%ai0, ai1, ai2 - Parameter aii and its T-derivatives (1xNC)
%Outputs:
%Parameters of the mixture, temperature and compositional derivatives

%Calculate the B mixture parameter
Bmix=zfeed*bii';

%Calculate the A mixture parameter
AX0=zfeed.*ai0;
AH0=AX0*kij;
AD1=ai0.*AH0;
Amix=zfeed*AD1';

%Calculate dA/dX
Amix_dX=2*AD1;

%Calculate derivatives if required
AX1=zfeed.*ai1;
AH1=AX1*kij;
Amix_dXdT=2*(ai0.*AH1+ai1.*AH0);
Amix_dT=0.5*(Amix_dXdT*zfeed');
dummydt=(ai2.*AH0+ai1.*AH1);
Amix_dT2=2*(dummydt*zfeed');


