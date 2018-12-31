function [p1,p2,origin] = findWorldFrameOrigin(ray1,ray2,O,back_plane)
%FINDWORLDFRAMEORIGIN Find the origin of the world frame, considering the
%midpoint of the intersection points between 2 rays of symmetric points and
%the plane z=back_plane (i.e. the plane containing the car plate)
%
    syms x1, syms y1, syms z1;
    syms x2, syms y2, syms z2;
    
    x0 = O(1); y0 = O(2); z0 = O(3);
    l1 = ray1(1); m1 = ray1(2); n1 = ray1(3);
    l2 = ray2(1); m2 = ray2(2); n2 = ray2(3);
    
    eqns = [ (x1-x0)*m1 - (y1-y0)*l1;
             (z1-z0)*m1 - (y1-y0)*n1;
             (x2-x0)*m2 - (y2-y0)*l2;
             (z2-z0)*m2 - (y2-y0)*n2;
             z1 - back_plane;
             z2 - back_plane;
           ];
    res = solve(eqns);
    
    p1(1) = double(res.x1);
    p1(2) = double(res.y1);
    p1(3) = double(res.z1);
    
    p2(1) = double(res.x2);
    p2(2) = double(res.y2);
    p2(3) = double(res.z2);
    
    origin = (p1+p2)/2;
end
