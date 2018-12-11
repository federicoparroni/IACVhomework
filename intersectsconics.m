function [r] = intersectsconics(coeffs1,coeffs2)
%intersectsconics Return the 4 intersection point between 2 conics 
%   returns [[x1,y1],[x2,y2],...] matrix
    A1 = coeffs1(1);
    B1 = coeffs1(2);
    C1 = coeffs1(3);
    D1 = coeffs1(4);
    E1 = coeffs1(5);
    F1 = coeffs1(6);
    
    A2 = coeffs2(1);
    B2 = coeffs2(2);
    C2 = coeffs2(3);
    D2 = coeffs2(4);
    E2 = coeffs2(5);
    F2 = coeffs2(6);
    
    syms x, syms y;
    eq1 = A1*x^2 + B1*x*y + C1*y^2 + D1*x + E1*y + F1;
    eq2 = A2*x^2 + B2*x*y + C2*y^2 + D2*x + E2*y + F2;
    eqns = [eq1, eq2];
    sol = solve(eqns);

    r = [double(sol.x), double(sol.y)];
end

