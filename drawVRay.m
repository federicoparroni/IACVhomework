function [] = drawVRay(M, O, p, style)
%DRAWVRAY Draw the viewing ray from viewpoint O and projecting into p of
%camera with rotation matrix M=K*R
%
    ray_direction = M \ p;
    
    x = [-1.5,1.5];
    ray = [ x;
        ray_direction(2)*(x-O(1))/ray_direction(1) + O(2);
        ray_direction(3)*(x-O(1))/ray_direction(1) + O(3)
      ];
    %plotLine3([ray(1,1),ray(3,1),-ray(2,1)], [ray(1,2),ray(3,2),-ray(2,2)],style,1);
    plotLine3([ray(1,1),ray(2,1),ray(3,1)], [ray(1,2),ray(2,2),ray(3,2)],style,1);
end

