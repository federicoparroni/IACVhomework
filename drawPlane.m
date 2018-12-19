function [] = drawPlane(coeff, color)
%DRAWPLANE Draw a 3d plane
%   
    a = coeff(1);
    b = coeff(2);
    c = coeff(3);
    d = coeff(4);

    [x,y] = meshgrid(-1:0.1:1);
    if c==0
        z = repelem(d,size(x,1)); %zeros(size(x, 1));
    else
        z = -1/c*(a*x + b*y + d);
    end
    
    surf(x,y,z, 'FaceColor',color);
end

