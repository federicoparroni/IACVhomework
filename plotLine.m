function [] = plotLine(a,b,c,support,style,width)
% plotLine Plot a line in the current figure
%   
    hold on;
    m = -a/b;
    q = -c/b;
    if style == 0
        plot(support, m*support + q, 'HandleVisibility','off', 'LineWidth',width);
    else
        plot(support, m*support + q, style, 'LineWidth',width);
    end
    hold off;
end

