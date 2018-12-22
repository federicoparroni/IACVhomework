function [] = plotLine(line,support,style,width)
% plotLine Plot a line in the current figure
% line: array representing line coefficients [a,b,c]
% support: x coords to draw
% style: plot style (see 'help plot'), e.g.:  'r-'
% width: width of the line
%   
    hold on;
    m = -line(1)/line(2);
    q = -line(3)/line(2);
    if style == 0
        plot(support, m*support + q, 'HandleVisibility','off', 'LineWidth',width);
    else
        plot(support, m*support + q, style, 'LineWidth',width);
    end
    hold off;
end

