function [] = plotPoint3(point, style)
%PLOTPOINT3 Draw a 3D point using the standard right hand reference
%   
    %plot3(-point(1),point(3),-point(2),style);
    plot3(point(1),point(2),point(3),style);
end

