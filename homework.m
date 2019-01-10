%%  IMAGE ANALYSIS & COMPUTER VISION HOMEWORK
%   Federico Parroni 10636497
%   

close all; clear all; clc;
warning('off', 'images:initSize:adjustingMag');

%% 1. Ellipsis detection
% Load the original image and compute the normalization factor and matrix
% to bring coordinates in the interval [0,1]
original_im = imread('Image1.jpg');

norm_f = max((size(original_im)));
norm_mx = diag([1/norm_f, 1/norm_f, 1]);

maxw = length(original_im) * 1.2;   % used to draw lines over the image

figure(); imshow(original_im); title('Original image');

%% 1.1 Preprocess the original image
% As first, enhance the contrast of the image by remapping pixels
% intensities. We will apply 2 different filters to detect both wheels of
% the car
gray_im = rgb2gray(original_im);
scale_factor = 4;
resized_gray_im = imresize(gray_im, 1/scale_factor);

curves = uint8([zeros(1,60) linspace(0,255,16) repelem(255,180)]);
figure(); plot(0:255, curves);
xlabel('original pixel intensity'); ylabel('transformed pixel intensity');
title('Preprocessing');

preprocessed_im = intlut(resized_gray_im, curves);
figure(); imshow(preprocessed_im); title('Preprocessed image');
%%
% Filter 1: binarize the preprocessed image with a global
% threshold computed using Otsu's method, then apply the 'remove' filter
% (set a pixel to 0 if its 4 neighbors are all 1, thus leaving only
% boundary pixels)
imf1 = imbinarize(preprocessed_im);
imf1 = bwmorph(imf1,'remove');
%%
% Filter 2: binarize again the original image, but with different
% threshold values computed by changing the neighborhood size
threshold1 = adaptthresh(preprocessed_im, 0.35, 'ForegroundPolarity','dark');
threshold2 = adaptthresh(preprocessed_im, 0.35, 'NeighborhoodSize',13);
thr1 = uint8(imbinarize(preprocessed_im, threshold1));
thr2 = uint8(imbinarize(preprocessed_im, threshold2));
%%
% Sum the thresholded images, in this way we can select white parts in
% all the thresholded images
imf2 = preprocessed_im + thr1 + thr2;
clear thr1; clear thr2;

% Dilate the filtered image to reduce black holes inside the wheel
imf2 = imbinarize(imf2);
imf2 = imdilate(imf2, strel('disk',2,8));

% Show the 2 filtered images
figure(); imshow(imf1); title('Filtered image 1');
figure(); imshow(imf2); title('Filtered image 2');

%% 1.2 Detect the ellipsis in the preprocessed images
% Get the ellipsis in both the filtered images using the regionprops function
s1 = regionprops(imf1, {'Centroid','MajorAxisLength','MinorAxisLength','Orientation'});
s2 = regionprops(imf2, {'Centroid','MajorAxisLength','MinorAxisLength','Orientation'});

% Filter all the ellipses detected before to exclude bad candidates
s1 = s1(arrayfun(@(x) x.MajorAxisLength < 120 && x.MajorAxisLength > 20 ...
                   && x.MinorAxisLength < 80 && x.MajorAxisLength > 20 ...
                   && (x.Orientation < -80 && x.Orientation > -90 ...
                   || x.Orientation < 78 && x.Orientation > 72), s1));
s2 = s2(arrayfun(@(x) x.MajorAxisLength < 120 && x.MajorAxisLength > 110, s2));

% Draw the ellipsis of interest
t = linspace(0,2*pi,60);
figure(); imshow(imf1); title('Ellipsis in filtered image 1');
for k=1:length(s1)
    drawEllipse(s1(k),t);
end
figure(); imshow(imf2); title('Ellipsis in filtered image 2');
for k=1:length(s2)
    drawEllipse(s2(k),t);
end

ellipsis = [s1; s2];
figure(); imshow(resized_gray_im); title('Ellipsis of interest');
for k=1:length(ellipsis)
    drawEllipse(ellipsis(k),t);
    % normalize the ellipsis
    ellipsis(k).Centroid = ellipsis(k).Centroid / norm_f * scale_factor;
    ellipsis(k).MajorAxisLength = ellipsis(k).MajorAxisLength / norm_f * scale_factor;
    ellipsis(k).MinorAxisLength = ellipsis(k).MinorAxisLength / norm_f * scale_factor;
end
clear s1; clear s2; clear t;

%% 2.1 Determine the ratio between the diameter of the wheel circles and the wheel-to-wheel distance 
% 1. Get the conic matrices associated to the ellipsis found before. We can
% retrieve the equation of the ellipsis from their centers, axis lengths
% and orientation.
%%
% 2. Find the tangent lines by intersecting the dual conics.
%%
% 3. Find two vanishing points from two pairs of parallel lines.
%%
% 4. Find the image of the line at infinity.
%%
% 5. Find the images of the circular points by intersecting ellipsis and
% line at infinity.
%%
% 6. Find the image of the degenerate dual conic by multiplying the images
% ofthe circular points.
%%
% 7. Find the rectification homography using SVD.

%% 2.1.1 Find the dual conics starting from the ellipsis
% The equation of an ellipse:
%%
%
% $$ ax^2 + by^2 + bxy + dx + ey + f = 0 $$
%
% can be expressed from its properties. If we know that the ellipse is
% center in $(Xc,Yc)$, with axis $a$ and $b$ and orientation $\theta$,
% where:
%
% $$ a = a^2 \sin^2 \theta + b^2 \cos^2 \theta $$
%
% $$ b = 2(b^2 - a^2)\sin \theta \cos \theta $$
%
% $$ c = a^2 \cos^2 \theta + b^2 \sin^2 \theta $$
%
% $$ d = -2ax_c -by_c $$
%
% $$ e = -bx_c -2cy_c $$
%
% $$ f = ax_c^2 +bx_cy_c +cy_c^2 -a^2b^2 $$

wheel_ant = ellipsis(2);
wheel_post = ellipsis(3);
wheel_ant_coeffs = ellipseCoeff(wheel_ant);
wheel_post_coeffs = ellipseCoeff(wheel_post);
% Find the dual conics, inverting the matrices of the conics
dual_wheel_ant_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(wheel_ant_coeffs)));
dual_wheel_post_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(wheel_post_coeffs)));

% Get the 4 tangent lines and plot them with different colors
tangent_lines = intersectsconics(dual_wheel_ant_coeffs,dual_wheel_post_coeffs);
% Plot the 4 tangent lines
figure(); imshow(gray_im);
styles=['r-','m-','b-','y-']; widths=[1 2 2 1];
for k=1:length(tangent_lines)
    line = [tangent_lines(k,:), 1] * norm_mx;
    plotLine(line,[0 maxw], styles(k), widths(k));
end
clear line; clear styles; clear widths; clear k;

%% 2.1.1 Intersect lines and conics
% Intersects lines and conics
uppertangentline = [tangent_lines(2,:), 1];
lowertangentline = [tangent_lines(3,:), 1];
upperwheelpoint1 = real(intersectsconics(wheel_post_coeffs, [0 0 0 uppertangentline]));
lowerwheelpoint1 = real(intersectsconics(wheel_post_coeffs, [0 0 0 lowertangentline]));
upperwheelpoint2 = real(intersectsconics(wheel_ant_coeffs, [0 0 0 uppertangentline]));
lowerwheelpoint2 = real(intersectsconics(wheel_ant_coeffs, [0 0 0 lowertangentline]));
% Take only one of the 2 identically tangent points (molteplicity 2)
upperwheelpoint1 = upperwheelpoint1(1,:);
lowerwheelpoint1 = lowerwheelpoint1(1,:);
upperwheelpoint2 = upperwheelpoint2(1,:);
lowerwheelpoint2 = lowerwheelpoint2(1,:);
% Get vertical lines through tangency points of each wheel
line_vert_post_wheel = cross([lowerwheelpoint1, 1],[upperwheelpoint1, 1]);
line_vert_ant_wheel = cross([lowerwheelpoint2, 1],[upperwheelpoint2, 1]);
% Get the vanishing points
vanish_pointz = cross(uppertangentline, lowertangentline);
vanish_pointz = vanish_pointz / vanish_pointz(3);
vanish_pointz = vanish_pointz.';
vanish_pointy = cross(line_vert_post_wheel, line_vert_ant_wheel);
vanish_pointy = vanish_pointy / vanish_pointy(3);
vanish_pointy = vanish_pointy.';
% Get the 4 intersection points
middle_point = cross( [tangent_lines(1,:), 1], [tangent_lines(4,:), 1] );
middle_point = middle_point / middle_point(3);
middle_line = cross( vanish_pointz, middle_point' ).';
post_wheel_points = intersectsconics(wheel_post_coeffs, [0 0 0 middle_line]);
ant_wheel_points = intersectsconics(wheel_ant_coeffs, [0 0 0 middle_line]);

% Plot all tangency and vanishing points and lines
hold on;
plot(upperwheelpoint1(1)*norm_f,upperwheelpoint1(2)*norm_f,'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint1(1)*norm_f,lowerwheelpoint1(2)*norm_f,'m^','MarkerSize',6,'HandleVisibility','off');
plot(upperwheelpoint2(1)*norm_f,upperwheelpoint2(2)*norm_f,'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint2(1)*norm_f,lowerwheelpoint2(2)*norm_f,'m^','MarkerSize',6,'HandleVisibility','off');
plot(middle_point(1)*norm_f,middle_point(2)*norm_f,'r^','MarkerSize',6,'HandleVisibility','off');
plot(vanish_pointz(1)*norm_f,vanish_pointz(2)*norm_f,'b*','MarkerSize',6,'HandleVisibility','off');
plot(vanish_pointy(1)*norm_f,vanish_pointy(2)*norm_f,'b*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(1,1)*norm_f,ant_wheel_points(1,2)*norm_f,'r*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(2,1)*norm_f,ant_wheel_points(2,2)*norm_f,'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(1,1)*norm_f,post_wheel_points(1,2)*norm_f,'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(2,1)*norm_f,post_wheel_points(2,2)*norm_f,'r*','MarkerSize',6,'HandleVisibility','off');
plotLine(line_vert_post_wheel * norm_mx, [0 maxw],'-',2);
plotLine(line_vert_ant_wheel * norm_mx, [0 maxw],'-',2);

% Get the image of the line at the infinity
inf_point1 = cross( lowertangentline, uppertangentline );
inf_point2 = cross( line_vert_post_wheel, line_vert_ant_wheel );
line_infinity = cross( inf_point1, inf_point2 );
plotLine(line_infinity * norm_mx, [0 maxw], '-', 2);
% Get the image of the circular points by intersecting the ellipsis and the
% line at the infinity
circ_point1 = intersectsconics(wheel_post_coeffs, [0 0 0 line_infinity]);
circ_point2 = intersectsconics(wheel_ant_coeffs, [0 0 0 line_infinity]);
circ_points = (circ_point1 + circ_point2) ./ 2;
legend('line1','line2','line3','line4', 'vertical1','vertical2', 'infinity');
hold off;

%% 2.1.2 Rectify the image using the image of the degenerate dual conic
% The image of the degenerate dual conic can be computed by multiplying the
% images of the circular points:
% 
% $\omega_{\infty} = I*J' + J*I'$
% 
I = [circ_points(1,:).'; 1];
J = [circ_points(2,:).'; 1];
image_dual_conic = I*J.' + J*I.';
% When we have w, we can find the homography using SVD
[~, DC, H] = svd(image_dual_conic);
normalization = sqrt(DC); normalization(3,3)=1;
H_n = normalization * H;
H = inv(H_n);

% Find the rectified intersection wheel points
a = H \ [post_wheel_points(2,:), 1]'; a = a/a(3);
b = H \ [post_wheel_points(1,:), 1]'; b = b/b(3);
c = H \ [ant_wheel_points(2,:), 1]'; c = c/c(3);
d = H \ [ant_wheel_points(1,:), 1]'; d = d/d(3);

% Compute the ratio distance - diameter from the rectified points
ratio = norm( a-b, 2) / norm( a-c, 2)

figure(); hold on;
plot(a(1),a(2),'or'); plot(b(1),b(2),'or');
plot(c(1),c(2),'ob'); plot(d(1),d(2),'ob');
hold off; title('Rectification of 4 wheel keypoints');


%% 2.2 Camera calibration
% We can compute K through the image of the absolute conic.
%%
% We can compute R through the mapping of the 3 orthogonal directions into
% the image vanishing points.

%% 2.2.1 Detect the features
% We first need the missing vanishing point.
% Get some features in the image using Harris method to find another pair
% of parallel lines.
figure(); imshow(original_im);
[feat1_x,feat1_y] = harris(gray_im, 4, 20);
hold on; plot(feat1_y,feat1_x,'r+'); hold off;

%%
% Take 4 feature points beloging to parallel lines in order to find the
% last vanishing point
stop_point1 =  [1549, 464, 1];
stop_point2 =  [1281, 440, 1];
plate_point1 = [1207, 1569, 1];
plate_point2 = [930, 1478, 1];

% Find the parallel lines through these pairs of points and normalize them
line_lights = cross( stop_point1 * norm_mx, stop_point2 * norm_mx );
line_plate = cross( plate_point1 * norm_mx, plate_point2 * norm_mx );

plotLine(line_lights * norm_mx, [-maxw maxw],'g-',2);
plotLine(line_plate * norm_mx, [-maxw maxw],'g-',2);
hold on;
vanish_pointx = cross( line_lights, line_plate );
vanish_pointx = vanish_pointx.' / vanish_pointx(3);
plot(vanish_pointx(1)*norm_f,vanish_pointx(2)*norm_f,'b*','MarkerSize',6,'HandleVisibility','off');

%% 2.2.2 Compute K from w
% Get the image of the absolute conic (w) from the constraints on
% orthogonal directions
syms fx, syms fy, syms u0, syms v0;
K = [fx, 0, u0; ...
     0, fy, v0; ...
     0, 0, 1];
w_dual = K*K.';
w = inv(w_dual);
eqs = [ vanish_pointz.' * w * vanish_pointy;
        vanish_pointz.' * w * vanish_pointx;
        vanish_pointy.' * w * vanish_pointx;
        [circ_points(1,:) 1] * w * [circ_points(1,:) 1].';
      ];
res = solve(eqs);
fx = real(double(res.fx)); fx = fx(fx > 0); fx = fx(1);
fy = real(double(res.fy)); fy = fy(fy > 0); fy = fy(1);
u0 = real(double(res.u0)); u0 = u0(1);
v0 = real(double(res.v0)); v0 = v0(1);
clear K, clear eqs, clear res;
% Print the matrix K (normalized)
K = [fx, 0, u0; ...
     0, fy, v0; ...
     0, 0, 1]


%% 2.3 3D Reconstruction
% We want to reconstruct the 3D location of some symmetric point in the
% image. We need the intrinsic and extrinsic parameters to do this.
% Find the location of some image keypoints in the scene using
% backprojection. We fix the world reference frame at the center of the car
% plate, with the y-axis contained in the vertical symmetry plane.
% We can find the rotation matrix of the camera with respect to the world
% reference frame. We can exploit the vanishing points to find the
% columns of R.

%%
% N.B.: in all the following plots, axis x, y, z are colored with red,
% green and blue respectively

%% 2.3.1 Fix the origin of reference frame in the camera frame origin
%

% Compute the rows of R and normalize each one:
r1 = K \ vanish_pointx; r1 = r1/norm(r1);
r3 = K \ vanish_pointz; r3 = r3/norm(r3);
% Last column can be found by cross product by the other two
r2 = cross(r3,r1);
plate_center = transpose([1033, 1513, 1] * norm_mx);

% Initially we set t = 0
t_0 = [0; 0; 0];

R = [r1 r2 r3]
P = K*[R t_0];
M = P(:,1:3);

O = null(P); O = O/O(4);

% Get some symmetric features
plate_left = transpose([901, 1465, 1] * norm_mx);
plate_right = transpose([1207, 1570, 1] * norm_mx);
% Use stop points found before
stop_left = transpose(stop_point1 * norm_mx);
stop_right = transpose(stop_point2 * norm_mx);
light_left = transpose([777, 936, 1] * norm_mx);
light_right = transpose([1521, 1116, 1] * norm_mx);

figure(); imshow(gray_im); hold on;
plot([901,1207,777,1521,stop_point1(1),stop_point2(1)], [1465,1570,936,1116,stop_point1(2),stop_point2(2)], 'r+');
hold off;

%%
% Find the backprojection rays of some symmetric points:
viewing_ray_platel = M \ plate_left;
viewing_ray_plater = M \ plate_right;
viewing_ray_lightl = M \ light_left;
viewing_ray_lightr = M \ light_right;
viewing_ray_stopl = M \ stop_left;
viewing_ray_stopr = M \ stop_right;

% Find the origin of the new reference frame
z_dist = 1;
[symm_point_platel,symm_point_plater,neworigin] = findWorldFrameOrigin(viewing_ray_platel,viewing_ray_plater,O, z_dist);
t_2 = neworigin.';

%%
% Knowning the viewing rays of the pairs of the projected
% image points that were symmetric in the scene, we can find the original
% 3D location by imposing that these points must be symmetric with respect
% to the symmetry plane found before. Everything is expressed in the
% reference frame W0 (camera frame rotated by R).
%
figure(); hold on;
plotPoint3([0,0,0],'wo',8);  % origin
% versors
plotLine3([0,0,0],[0.5,0,0],'r',2);
plotLine3([0,0,0],[0,0.5,0],'g',2);
plotLine3([0,0,0],[0,0,0.5],'b',2);
% viewpoint
plotPoint3(O,'k^',8);
% new origin in the symmetry plane
plotPoint3(neworigin,'r^',8);

drawVRay(M,O,plate_left, 'c', 2);
drawVRay(M,O,plate_right, 'c', 2);

drawVRay(M,O,light_left, 'm', 2);
drawVRay(M,O,light_right, 'm', 2);

drawVRay(M,O,stop_left, 'g', 2);
drawVRay(M,O,stop_right, 'g', 2);

% planes z=1 and symmetry plane
drawPlane([0; 0; 1; -z_dist], 'k', 3);
drawPlane([1; 0; 0; -neworigin(1)], 'y', 2);

% Solve the system to find the pairs of symmetric features
[symm_point_lightl,symm_point_lightr] = findSymmetricFeatures(viewing_ray_lightl, viewing_ray_lightr, O, neworigin(1));
[symm_point_stopl,symm_point_stopr] = findSymmetricFeatures(viewing_ray_stopl, viewing_ray_stopr, O, neworigin(1));

% Plot the 3D locations of the symmetric features
plotPoint3(symm_point_platel,'cp',8);
plotPoint3(symm_point_plater,'cp',8);
plotPoint3(symm_point_lightl,'m^',6);
plotPoint3(symm_point_lightr,'m^',6);
plotPoint3(symm_point_stopl,'gs',6);
plotPoint3(symm_point_stopr,'gs',6);

plotCamera('Location',O(1:3),'Orientation',R,'Size',0.1);

pbaspect([1,1,1]);
grid on; title('3D Reconstruction (in the reference frame W0)');
xlabel('x'); ylabel('y'); zlabel('z');
legend('origin','x axis','y axis','z axis', 'viewpoint', 'car frame origin', ...
         'left plate ray','right plate ray','left light ray', ...
         'right light ray','left stop ray','right stop ray', ...
         'car back plane','symmetry plane', ...
         'left plate','right plate','left light','right light', ...
         'left stop','right stop');
set(gca, 'CameraPosition', [0,1,-1]);
hold off;

%% 2.3.2 Fix the reference frame in the car symmetry plane
% We can express all the points in the reference frame W1 (with origin in
% the symmetry plane of the car). We need to build the matrix T1 (that is a
% pure translation of -neworigin).
T1 = [eye(3) -neworigin; 0 0 0 1];

symm_point_platel_W1 = T1*[symm_point_platel; 1];
symm_point_plater_W1 = T1*[symm_point_plater; 1];
symm_point_lightl_W1 = T1*[symm_point_lightl; 1];
symm_point_lightr_W1 = T1*[symm_point_lightr; 1];
symm_point_stopl_W1 = T1*[symm_point_stopl; 1];
symm_point_stopr_W1 = T1*[symm_point_stopr; 1];

% Replot all points in W1
figure();
hold on;
plotPoint3([0,0,0],'wo',8);
% versors
plotLine3([0,0,0],[0.1,0,0],'r',2);
plotLine3([0,0,0],[0,0.1,0],'g',2);
plotLine3([0,0,0],[0,0,0.1],'b',2);
% viewpoint and symmetric points
plotPoint3(symm_point_platel_W1,'cp',8);
plotPoint3(symm_point_plater_W1,'cp',8);
plotPoint3(symm_point_lightl_W1,'m^',6);
plotPoint3(symm_point_lightr_W1,'m^',6);
plotPoint3(symm_point_stopl_W1,'gs',6);
plotPoint3(symm_point_stopr_W1,'gs',6);

pbaspect([1,1,1]);
grid on; title('3D Reconstruction (in the reference frame W1)');
legend('origin','x axis','y axis','z axis', ...
         'left plate','right plate','left light','right light', ...
         'left stop','right stop');
set(gca, 'CameraPosition', [0,1,-1]);

%% 2.4 Localize the camera in the world frame
% Now that we have fix the reference frame W1 with respect to W0, we can
% find the location of the camera in this new frame.

camera_location = T1*O
plot3(camera_location(1), camera_location(2), camera_location(3), ...
   'k^', 'MarkerSize',8, 'DisplayName','camera');

