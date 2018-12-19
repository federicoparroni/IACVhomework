function [] = plotLine3(p1,p2, style, width)
%PLOTLINE3 Draw a 3D line starting from p1 and ending at p2
%   
	plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)], style, 'LineWidth',width);

end

