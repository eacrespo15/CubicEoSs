function[data]=get_module_database()
%Database file of critical properties and binary interaction parameters
%to use with the cEoSs
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%[data] is a structure with the following fields:
%name - a short name for each component
%TC - critical temperature (K)
%PC - critical pressure (MPa)
%ACENTRIC - acentric factor
%binary - Matrix of binary interaction parameters (NCxNC), where NC is the
%number of components in the database

%Name of the components
data.name=...
    {'C1','C2','C3','IC4','C4'...
    'IC5','C5','C6','C7','C8'...
    'C9','H2O','N2','CO2','H2S'};
%Critical Temperature (K)
data.TC=[...
190.6 305.4 369.8 408.1 425.2...
460.4 469.6 507.4 540.2 568.8...
594.6 647.3 126.2 304.2 373.2];
%Critical Pressure (MPa)
data.PC=[...
    4.599 4.8827 4.2445 3.6468 3.7988...
    3.3834 3.3733 2.9681 2.7351 2.4819...
    2.3096 22.0834 3.3936 7.3746 8.9374];
%Acentric Factor
data.ACENTRIC=[...
    0.008 0.098 0.152 0.176 0.193...
    0.227 0.251 0.296 0.351 0.394...
    0.440 0.344 0.040 0.225 0.100];
%Binary Parameters
NCD=length(data.PC);
data.binary=zeros(NCD,NCD);
data.binary(1,12:15)=[0.45 0.02 0.12 0.08];
data.binary(2,12:15)=[0.45 0.06 0.15 0.07];
data.binary(3,12:15)=[0.53 0.08 0.15 0.07];
data.binary(4,12:15)=[0.52 0.08 0.15 0.06];
data.binary(5,12:15)=[0.52 0.08 0.15 0.06];
data.binary(6,12:15)=[0.50 0.08 0.15 0.06];
data.binary(7,12:15)=[0.50 0.08 0.15 0.06];
data.binary(8,12:15)=[0.50 0.08 0.15 0.05];
data.binary(9,12:15)=[0.50 0.08 0.15 0.04];
data.binary(10,12:15)=[0.50 0.08 0.15 0.04];
data.binary(11,12:15)=[0.50 0.08 0.15 0.03];
data.binary(14,15)=0.12;
for i=1:NCD-1
    for j=1:NCD
        data.binary(j,i)=data.binary(i,j);
    end
end

end

