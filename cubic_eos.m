function [phasefound,Z,fug,fugT,fugP,fugX,AUX] = cubic_eos(NC,T,P,zfeed,phasetype,solvertype,SQTC,kij,delta1,delta2,bci,aci,mfunc)
%Core function for the calculations using cEoSs
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Inputs:
%NC         - Number of components
%indexes    - Indexes of the system components in the database
%T          - Temperature (K)
%P          - Pressure(atm)
%zfeed      - Mixture's composition
%EoS        - Type of Cubic EoS (0-SRK || 1- Peng-Robinson)
%phasetype  - 1:Liquid || -1:Vapor || 0-Min. Gibbs Energy
%solvertype - 1:Mollerup || 2- Mollerups with roots sorting
%
%Outputs:
%phasefound- phase found by the EoS: If phase type= 1/-1 it returns 1 if the
%   requested phase was found and -1 otherwise. If phasetype=0 it returns 2 if a
%   liquid root was found or -2 if a vapor root was found!
%Zcalc - Compressibility factor 
%fug - Vector of ln(fugacity coefficients) (1xNC)
%fugT - T-derivatives of fug (1xNC)
%fugP - P-derivatives of fug (1xNC)
%fugX - Scaled composition derivatives of fug (NCxNC)
%AUX - Various residual properties Hres, Sres, Cpres,...


%Initialize the variable containing various properties
AUX=zeros(1,10);

%Get T-dependent terms of the EoS
[Ai0,Ai1,Ai2]=set_temperature(T,SQTC,aci,mfunc,NC);

%Get mixture parameters A and B and their appropriate derivatives
[Bmix,Amix,Amix_dT,Amix_dT2,Amix_dX,Amix_dXdT]=mixparam(zfeed,kij,bci,Ai0,Ai1,Ai2);

%Solve the Cubic EoS in terms of Compressibility factor
TINV=1/T;
POVERT=P*TINV;
APT=Amix*POVERT;
BPT=Bmix*POVERT;

[Z,ZV2,dG]=cubic_solver(phasetype,APT,BPT,delta1,delta2,solvertype);


%Calculate the Volume
PINV=1/P;
V=Z*T*PINV; %This is V/R 

%It is worth to remember that the whole EoS was solved eliminating the gas constant

%Return the appropriate phase found indicator
phasefound=1;
if (V>=3*Bmix)
    phasefound=-1;
end

if phasetype ~=0
    phasefound=phasefound*phasetype;
else
    phasefound=2*phasetype;
end

%Auxiliary dummy values
BINV=1/Bmix;
S1=1/(V+delta1*Bmix);
S2=1/(V+delta2*Bmix);

%Calculation of partial pressure derivatives
%P1=P*TINV+Amix*S1*S2;          %This is equivalent to 1/(V-b) of eq. 5
P1=1/(V-Bmix);                  %First term of eq.5
PA=-S1*S2;                      %d(P/T)/d(A/T)
P2=Amix*PA;                     %Equivalent to the second term of eq.5
dummy1=delta1*S1+delta2*S2;     %Useful value 
PB=P1*P1-dummy1*P2;             %This is the derivative of P/T in order to Bmix
PN=P1;

%Calculation of Helmholtz derivatives and useful dummy values
F1=log(V*P1);                       %(FN) eq. 75 Chapter 3
dummy2=log(S1/S2)/(delta2-delta1);  %(XL2)Useful value
f=dummy2*BINV;                      %(-FA)f from eq. 62
F2=Amix*f;                          %Second term from eq. 60
FF=F1-F2;                           %F function from eq. 60 - Helmholtz energy/RT
dummy3=-V*PA;                       %(GB) This is equal to 1/(delta1-delta2) * df2/db onde f2=ln(V+delta1*B)/(V+Delta2*B)
F1B=P1;                             %(FnB) Derivative of F1 in order to B;
F2B=(Amix*dummy3-F2)*BINV;               %Derivative of F2 in order to B;
FB=F1B-F2B;                         %Derivative of F in order to B
fB=F2B/Amix;                        %(-FAB) Derivative of f in order to B
dummy4=-dummy3*dummy1;              %(GBB) Necessary to calculate FBB
F2BB=(Amix*dummy4-2*F2B)*BINV;      %Derivative of F2B in order to B;
F1BB=P1*P1;                         %eq.97
FBB=F1BB-F2BB;                      %Second-order derivative of F in order to B


%Calculation of residual Entropy, Enthalpy, and Gibbs Energy
dummy5=log(Z);
dFdT=-f*Amix_dT;
SR=-T*dFdT-FF+dummy5; %This is the value of Sres/R
HR=T*(-T*dFdT+Z-1);   %This is the value of Hres/R
GR=HR-T*SR;           %This is the value of Gres/R

%Calculation of useful derivatives of F or P
DPDV=-P1*P1-P2*(S1+S2); %dP/dV (or more precisely d(P/T)/d(V/R))
DPDT=T*PA*Amix_dT+P/T;      %Note that d2F_dTdV=-PA*Amix_dT
DVDT=-DPDT/DPDV*TINV;   %d(V/R)/dT
d2FdT2=-f*Amix_dT2;     

CVR=-T*(T*d2FdT2+2*dFdT); %Cv/R
CPR=CVR-DPDT*DPDT/DPDV-1; %Cp/R
dHdP=V-T*DVDT;            %dH/dP
dSdP=-DVDT+PINV;          %Why is 1/P added here?


%Calculation of ln(fugcoef)
FNP=F1-dummy5;
fug=FNP+FB*bci-f*Amix_dX;

%Calculation of T and P derivatives of ln(fugcoef)
dPdX=PN+PB*bci+PA*Amix_dX; %The dPdX value is divided by T
dVdX=-dPdX/DPDV;           %Because DPDV is d(P/T)... This is Vi/R
fugP=dVdX*TINV-PINV;       
d2FdTdX=-fB*bci*Amix_dT-f*Amix_dXdT;
fugT=d2FdTdX+TINV-DPDT*TINV*dVdX;

%Calculation of Compositional derivatives of ln(fugcoef)
fugX=zeros(NC,NC);
for i=1:NC
    for j=1:NC
        T1=F1B*(bci(i)+bci(j));
        T2=-fB*(bci(i)*Amix_dX(j)+bci(j)*Amix_dX(i));
        T3=0; %If the expression for B is modified this has to be modified
        T4=FBB*bci(i)*bci(j);
        T5=-2*f*Ai0(i)*Ai0(j)*kij(i,j);
        d2FdX2=T1+T2+T3+T4+T5;
        T6=dPdX(i)*dPdX(j)/DPDV;
        fugX(i,j)=d2FdX2+1+T6;
    end
end


%Store auxiliary properties
AUX(1)=ZV2;
AUX(2)=dG;
AUX(3)=HR;
AUX(4)=SR;
AUX(5)=CPR;
AUX(6)=dSdP;
AUX(7)=dHdP;
AUX(8)=DPDV*T;
AUX(9)=DVDT;
AUX(10)=DPDV;
end

