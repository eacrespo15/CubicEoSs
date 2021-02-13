function [AC0,AC1,AC2] = set_temperature(T,ALF,AC,m,nc)
%This funcion calculates the temperature dependent part of the EoS
%pure-component parameters, namely sqrt(A/RT) but as we eliminated the R
%constant AC0 represents sqrt(~A/T).
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Inputs:
%T          - Temperature (K)
%ALF        - 1/sqrt(Tc)                        (1xNC)
%AC         - Vector of ac(i) values            (1*NC)
%m          - m function for each component     (1xNC)
%nc         - Number of components
%
%Outputs:
%AC0= sqrt(A/T)
%AC1= d(AC0)/dT
%AC2= d2(AC0)/dT2

%Initialize the output variables
AC0=zeros(1,nc);
AC1=zeros(1,nc);
AC2=zeros(1,nc);

%Useful dummy variables
SQTR=1/sqrt(T);
T2R=0.5/T;
T2RF=-3*T2R;

%Calculation of the output variables
for i=1:nc
    Q1=AC(i)*(1+m(i))*SQTR;
    AC0(i)=Q1-AC(i)*m(i)*ALF(i);
    AC1(i)=-Q1*T2R;
    AC2(i)=T2RF*AC1(i);
end
end

