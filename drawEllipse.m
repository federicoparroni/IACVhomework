function [] = drawEllipse(s, t)
% drawEllipse Draw an ellipse.
%   s: structure returned by regionprops
%   t: support from 0 to 2pi, like: linspace(0,2*pi,60)
    
    a = s.MajorAxisLength/2;
    b = s.MinorAxisLength/2;
    Xc = s.Centroid(1);
    Yc = s.Centroid(2);
    phi = deg2rad(-s.Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    
    hold on;
    plot(x,y,'g','Linewidth',2,'HandleVisibility','off');
    hold off;
    
end

