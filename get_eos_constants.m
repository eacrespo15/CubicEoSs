function [delta1,delta2,CA,CB] = get_eos_constants(c)
%This subroutines calculates omegaB and sqrt(omegaA) that are used
%to calculate the ac and bc constants of the SRK/PR EoSs
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Reference: M. L. Michelsen and J.M. Mollerup Thermodynamic Models:
%Fundamental Aspects and Computational Aspects Cap.3
%
%Notes:
%c=0 yields the Redlich-Kwong EoS
%c=1 yields the Peng-Robinson EoS

%Calculation of delta 1 and delta 2
delta1=0.5*(c+1+sqrt((c+1)*(c+1)+4*c)); %Eq.10 pag.75
delta2=-c/delta1;                       %Eq. 9 pag.75
d1=1+c;                                 %d1=delta1+delta2=1+c
d2=-c;                                  %d2=delta1*delta2=-c

dummy1=power(1+delta1,1/3);
dummy2=power(1+delta2,1/3);
y=1+dummy1*dummy2*(dummy1+dummy2);      %Eqs. 16-18

dummy3=(3*y*(y+d1)+d1*d1-d2)/(3*y+c);   % Eq. 20

Zc=y/(3*y+c);                           %Eq. 19

OMEGAB=Zc/y;                            %Eq. 19 and Eq.24
OMEGAA=dummy3*OMEGAB;                   %Eq. 23

%For convenience the following variables are output instead of the omega
%values
CA=sqrt(OMEGAA);
CB=OMEGAB;
end

