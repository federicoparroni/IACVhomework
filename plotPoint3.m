function [] = plotPoint3(point, style, size)
%PLOTPOINT3 Convenience function to draw a 3D point with graphical style
%   
    plot3(point(1),point(2),point(3),style, 'MarkerSize',size);
end

