function [SQTC,kij,delta1,delta2,bci,aci,mfunc] = init_cubic(NC,EQT,INDEX)
%Initialization routine for the Cubic EoSs
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Inputs:
%NC    - Number of Components
%EQT   - Equation type (0-SRK, 1-PR)
%INDEX - Database indexes of the components present in the system
%
%Outputs:
%SQTC  - inverse of the square root of the critical temperatures (1xNC)
%kij   - Matrix of binary interaction parameters, already in the form (1-kij) (NC x NC)
%delta1 and delta 2- Constants for the EoS
%bci   - bc(i) values for every component (1xNC)
%aci   - sqrt(ac(i)) values for every component (1xNC)
%mfunc - Values of mSRK or mPR for every component (1xNC)

%Read the database
[data]=get_module_database();

%Get the critical properties
Tc=data.TC(INDEX);
Pc=data.PC(INDEX);
W=data.ACENTRIC(INDEX);

%Store useful variables
SQTC=1./sqrt(Tc);

%Get the matrix of binary interaction parameters in the form (1-kij)
[kij]=get_binary(data.binary,INDEX,NC);

%Read the EoS specific constants
if EQT==0
    C=0; %SRK
else
    C=1; %PR
end
[delta1,delta2,CA,CB]=get_eos_constants(C);

%Calculate ac and bc for each component (eqs. 23 and 24 of Cap.3)
aci=zeros(1,NC);
bci=zeros(1,NC);
mfunc=zeros(1,NC);
for i=1:NC
    bci(i)=CB*Tc(i)/Pc(i);
    aci(i)=CA*Tc(i)/sqrt(Pc(i));
    if C==0
        mfunc(i)=0.480+W(i)*(1.574-0.175*W(i));
    else
        mfunc(i)=0.37464+W(i)*(1.54226-0.26992*W(i));
    end
end

end

