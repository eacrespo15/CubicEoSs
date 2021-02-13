function[parameters] = get_binary(matrix,index,NC)
%This function retrives the binary interaction parameters matrix in the
%form (1-kij) for the system componentes.
%
%Emanuel A. Crespo
%PhD in Chemical Engineering: 2017-2021
%University of Aveiro/CICECO Aveiro Institute of Materials
%
%Code last revised in: February 2021
%
%Inputs:
%matrix- kij matrix from database
%index-  index of the different components present in the system
%NC- Number of components
%
%Outputs:
%parameters- matrix (NCxNC) containing the 1-kij values for all the systems
%sub-binary systems.

parameters=ones(NC,NC);
for i=1:NC
    for j=1:NC
        parameters(i,j)=1-matrix(index(i),index(j));
        parameters(j,i)=parameters(i,j);
    end
end
end

