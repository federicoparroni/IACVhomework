function [] = drawPlane(coeff, color, size)
%DRAWPLANE Draw a 3d plane
%   
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);

    if c==0
        [y,z] = meshgrid(-size:0.1:size);
        x = -1/a*(b*y + c*z + d);
    else
        [x,y] = meshgrid(-size:0.1:size);
        z = -1/c*(a*x + b*y + d);
    end
    
    surf(x,y,z, 'FaceColor',color);
end

