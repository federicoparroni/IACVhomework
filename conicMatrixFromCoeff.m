function [matrix] = conicMatrixFromCoeff(coeffs)
% conicMatrixFromCoeff Return the conic matrix from its coefficients
%   
    A = coeffs(1);
    B = coeffs(2);
    C = coeffs(3);
    D = coeffs(4);
    E = coeffs(5);
    F = coeffs(6);
    
    matrix(1,1) = A;
    matrix(1,2) = B/2;
    matrix(2,1) = B/2;
    matrix(1,3) = D/2;
    matrix(3,1) = D/2;
    matrix(2,2) = C;
    matrix(2,3) = E/2;
    matrix(3,2) = E/2;
    matrix(3,3) = F;
end

