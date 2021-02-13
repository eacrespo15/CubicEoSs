function [Z,ZV2,dG] = cubic_solver(phasetype,A,B,delta1,delta2,solvertype)
%This subroutine solves the cubic EoS for Z
%Searches for multiple roots
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%References: 1- https://www.e-education.psu.edu/png520/m11_p6.html
%2- Michelsens and Mollerups book
%3- Claudio Olivera-Fuentes, "The Optimal Solution of Cubic Equations of
%State", Latin American Applied Research, 23:243-256 (1993)
%
%Inputs:
%phasetype: 1:Liquid || -1: Vapor || 0: Min. Gibbs Energy
%A: AmixP/T
%B: BmixP/T
%delta1 and delta2: Parameters from the cubic EoS
%solvertype:    1-Mollerup's method
%               2-Mollerup's method with root order and selection

toldz=1E-7; %Tolerance for deltaZ
tolqz=1E-7; %Tolerance for the Q Function
if solvertype==1
    [Z,ZV2,dG]=cubic_solver1(phasetype,A,B,delta1,delta2,toldz,tolqz);
elseif solvertype==2
    [Z,ZV2,dG]=cubic_solver2(phasetype,A,B,delta1,delta2,toldz,tolqz);
else
    error('Invalid solver type selected for Cubic EoS');
end

end
    
function [Z,ZV2,dG] = cubic_solver1(phasetype,A,B,delta1,delta2,toldz,tolqz)
%This subroutine solves the cubic EoS for Z by the Mollerup's method
%Searches for multiple roots
%
%References: 1- https://www.e-education.psu.edu/png520/m11_p6.html
%2- Michelsens and Mollerups book
%3- Claudio Olivera-Fuentes, "The Optimal Solution of Cubic Equations of
%State", Latin American Applied Research, 23:243-256 (1993)
%
%Inputs:
%phasetype: 1:Liquid || -1: Vapor || 0: Min. Gibbs Energy
%A: AmixP/T
%B: BmixP/T
%delta1 and delta2: Parameters from the cubic EoS

%Initialize output variables
ZV2=0;
dG=0;

%Terms from the cubic polynomial in Z
dummy1=delta1+delta2;
dummy2=delta1*delta2;
dummy3=B*B;

d2=B*(dummy1-1)-1;
d1=A-B*dummy1+dummy3*(dummy2-dummy1);
d0=-(A*B+(dummy3+B*dummy3)*dummy2);

%Initial guess for Z
if phasetype==1
    Z=B; %which in reality is bP/RT
elseif phasetype==-1
    Z=1; %vapour-like estimate
else
    Z=-d2/3;
    Q=Z*((Z+d2)*Z+d1)+d0; %Q is the cubic EoS written in the form of a polynomial for Z
    Z=B;
    if Q<0
        Z=Z+1;
    end
end

%Setting up Third order Newton for Cubic EoS to converge first root
deltaZ=1;
Q=1;
while abs(deltaZ)>toldz && abs(Q)>tolqz
    DQ2=3*Z+d2;                         %Note that the q''(z) is divided by two here to simplify later
    DQ1=Z*(DQ2+d2)+d1;                  %q'(z) is calculated here using the q''/2
    Q=Z*((Z+d2)*Z+d1)+d0;               %q(z)
    IDQ1=1/DQ1;                         %1/q'(z)
    dummy4=Q*IDQ1;                      %q(z)/q'(z)
    deltaZ=dummy4*(1+DQ2*IDQ1*dummy4);  %This is delta but is lacking the minus sign
    Z=Z-deltaZ;                         %Because the minus was missing above here we do Z=Z-deltaZ
end

%Factorize the polynomial by knowing the first root x1
%So q(Z) may be rewritten as q(Z)=(Z-x1)(Z^2+F*Z+G)=0
F=Z+d2;
G=Z*F+d1;
D=F*F-4*G; %Expression inside the sqrt of quadratic expression

if D>=0 %If D < 0 there are no more real roots of the polynomial
    Z1=0.5*(abs(F)+sqrt(D));
    if F>0
        Z1=-Z1;
    end
    if Z>Z1
        Z1=G/Z1;
    end
    
    %Refinement of solution by a single Newton Step
    DQ2=3*Z1+d2;                            %Note that the q''(z) is divided by two here to simplify later
    DQ1=Z1*(DQ2+d2)+d1;                     %q'(z) is calculated here using the q''/2
    Q=Z1*((Z1+d2)*Z1+d1)+d0;                %q(z)
    IDQ1=1/DQ1;                             %1/q'(z)
    dummy5=Q*IDQ1;                          %q(z)/q'(z)
    deltaZ1=dummy5*(1+DQ2*IDQ1*dummy5);     %This is delta but is lacking the minus sign
    Z1=Z1-deltaZ1;                          %Because the minus was missing above here we do Z=Z-deltaZ
    
    if Z1>=B
        if phasetype==0
            %Calculate Excess Gibbs Energy Difference For Multiple Solutions
            d1B=delta1*B;
            d2B=delta2*B;
            F=log((Z-B)/(Z1-B))+A/(B*(delta2-delta1))*log((Z+d2B)*(Z1+d1B)/(Z+d1B)/(Z1+d2B));
            dG=abs(F);
            ZV2=Z1;
            if F<0       %Expression for Excess Gibbs Energy Difference (search for it)
                ZV2=Z;
                Z=Z1;
            end
        elseif phasetype==1 %If we want the liquid phase we store in Z the root with lower Z
            if Z1<Z
                Z=Z1;
            end
        else
            if Z1>Z
                Z=Z1;       %If we want the vapor phase we store in Z the root with higher Z
            end
        end
    end
end   
end

function [Z,ZV2,dG] = cubic_solver2(phasetype,A,B,delta1,delta2,toldz,tolqz)
%This subroutine solves the cubic EoS for Z by the Mollerup's method
%Searches for multiple roots
%
%References: 1- https://www.e-education.psu.edu/png520/m11_p6.html
%2- Michelsens and Mollerups book
%3- Claudio Olivera-Fuentes, "The Optimal Solution of Cubic Equations of
%State", Latin American Applied Research, 23:243-256 (1993)
%
%Inputs:
%phasetype: 1:Liquid || -1: Vapor || 0: Min. Gibbs Energy
%A: AmixP/T
%B: BmixP/T
%delta1 and delta2: Parameters from the cubic EoS

%Initialize output variables
ZV2=0;
dG=0;

%Terms from the cubic polynomial in Z
dummy1=delta1+delta2;
dummy2=delta1*delta2;
dummy3=B*B;

d2=B*(dummy1-1)-1;
d1=A-B*dummy1+dummy3*(dummy2-dummy1);
d0=-(A*B+(dummy3+B*dummy3)*dummy2);

%Initial guess for Z
if phasetype==1
    Z=B; %which in reality is bP/RT
elseif phasetype==-1
    Z=1; %vapour-like estimate
else
    Z=-d2/3;
    Q=Z*((Z+d2)*Z+d1)+d0; %Q is the cubic EoS written in the form of a polynomial for Z
    Z=B;
    if Q<0
        Z=Z+1;
    end
end

%Setting up Third order Newton for Cubic EoS to converge first root
deltaZ=1;
while abs(deltaZ)>toldz && abs(Q)>tolqz
    DQ2=3*Z+d2;             %Note that the q''(z) is divided by two here to simplify later
    DQ1=Z*(DQ2+d2)+d1;      %q'(z) is calculated here using the q''/2
    Q=Z*((Z+d2)*Z+d1)+d0;   %q(z)
    IDQ1=1/DQ1;             %1/q'(z)
    dummy4=Q*IDQ1;          %q(z)/q'(z)
    deltaZ=dummy4*(1+DQ2*IDQ1*dummy4); %This is delta but is lacking the minus sign
    Z=Z-deltaZ;             %Because the minus was missing above here we do Z=Z-deltaZ
end

%Factorize the polynomial by knowing the first root x1
%So q(Z) may be rewritten as q(Z)=(Z-x1)(Z^2+F*Z+G)=0
F=Z+d2;
G=Z*F+d1;
D=F*F-4*G; %Expression inside the sqrt of quadratic expression

%Calculate the other two roots if existing
if D>=0
    Z2=(-F+sqrt(D))/2;
    Z3=(-F-sqrt(D))/2;
    raizes=[Z,Z2,Z3];
    raizes=raizes(raizes>0);
    raizes=sort(raizes);
    if phasetype==0
        %Calculate Excess Gibbs Energy Difference For Multiple Solutions
        d1B=delta1*B;
        d2B=delta2*B;
        F=log((raizes(1)-B)/(raizes(end)-B))+A/(B*(delta2-delta1))*log((raizes(1)+d2B)*(raizes(end)+d1B)/(raizes(1)+d1B)/(raizes(end)+d2B));
        dG=abs(F);
        ZV2=raizes(end);
        if F<0       %Expression for Excess Gibbs Energy Difference (search for it)
            ZV2=Z;
            Z=raizes(end);
        end
    elseif phasetype==1 %If we want the liquid phase we store in Z the root with lower Z
        if raizes(end)<Z
            Z=raizes(end);
        end
    else
        if raizes(end)>Z
            Z=raizes(end);       %If we want the vapor phase we store in Z the root with higher Z
        end
    end
    
end
end
