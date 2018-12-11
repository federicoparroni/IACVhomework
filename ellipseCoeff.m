function [coefficients] = ellipseCoeff(ellipse)
% ellipseCoeff Returns the coefficients of the equation of the ellipse
% starting from its parameters
%   ellipse: struct returned by regionprops
%
%   Returns: [1]x^2 + [2]y^2 + [3]xy + [4]x + [5]y + [6]
    Xc = ellipse.Centroid(1);
    Yc = ellipse.Centroid(2);
    a = ellipse.MajorAxisLength/2;
    b = ellipse.MinorAxisLength/2;
    phi = deg2rad(-ellipse.Orientation);

    A = (a^2)*sin(phi)^2 + (b^2)*cos(phi)^2;
    B = 2*(b^2 - a^2)*sin(phi)*cos(phi);
    C = (a^2)*cos(phi)^2 + (b^2)*sin(phi)*sin(phi);
    D = -2*A*Xc -B*Yc;
    E = -B*Xc -2*C*Yc;
    F = A*(Xc^2) +B*Xc*Yc +C*(Yc^2) -(a^2)*(b^2);
    
    coefficients(1) = A;
    coefficients(2) = B;
    coefficients(3) = C;
    coefficients(4) = D;
    coefficients(5) = E;
    coefficients(6) = F;
    
end

