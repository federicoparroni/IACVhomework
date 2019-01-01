function [p1,p2] = findSymmetricFeatures(ray1,ray2,O, symm_plane)
%FINDSYMMETRICFEATURES Find 2 points belonging to the rays, passing through
% O, symmetric with respect to the plane x=symm_plane
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
             (x1 + x2)/2 - symm_plane;
             z1 - z2;
           ];
    res = solve(eqns);
    
    p1 = [double(res.x1);
          double(res.y1);
          double(res.z1) ];
    
    p2 = [double(res.x2);
          double(res.y2);
          double(res.z2) ];
end

