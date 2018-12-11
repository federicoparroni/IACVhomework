function [coeffs] = conicCoeffFromMatrix(matrix)
% conicCoeffFromMatrix Return the conic coefficients from its matrix
%   
    A = matrix(1,1);
    B = matrix(1,2) * 2;
    D = matrix(1,3) * 2;
    C = matrix(2,2);
    E = matrix(2,3) * 2;
    F = matrix(3,3);
    
    coeffs = [A B C D E F];
end

