%%  IMAGE ANALYSIS & COMPUTER VISION HOMEWORK
%   Federico Parroni 10636497
%   

close all; clear all; clc;

%% 1. Ellipsis detection
% Load the original image, resize for faster computations and convert it
% to grayscale
scale_factor = 4;
original_im = rgb2gray(imresize(imread('Image1.jpg'), 1/scale_factor));
%im_width = length(original_im);
norm_f = max((size(original_im)))*scale_factor;
norm_mx = diag([1/norm_factor, 1/norm_factor, 1]);

maxw = length(original_im) * 1.2;   % used to draw lines over the image

figure(); imshow(original_im); title('Original image');

%% 1.1 Preprocess the original image
% As first, enhance the contrast of the image by remapping pixels
% intensities. We will apply 2 different filters to detect both wheels of
% the car
curves = uint8([zeros(1,60) linspace(0,255,16) repelem(255,180)]);
figure(); plot(0:255, curves);
xlabel('original pixel intensity'); ylabel('transformed pixel intensity');
title('Preprocessing');

preprocessed_im = intlut(original_im, curves);
figure(); imshow(preprocessed_im); title('Preprocessed image');

% Filter 1: binarize the preprocessed image with a global
% threshold computed using Otsu's method, then apply the 'remove' filter
% (set a pixel to 0 if its 4 neighbors are all 1, thus leaving only
% boundary pixels)
imf1 = imbinarize(preprocessed_im);
imf1 = bwmorph(imf1,'remove');

% Filter 2: binarize again the original image, but with different
% threshold values computed by changing the neighborhood size
threshold1 = adaptthresh(original_im, 0.35, 'ForegroundPolarity','dark');
threshold2 = adaptthresh(original_im, 0.35, 'NeighborhoodSize',13);
t1 = uint8(imbinarize(original_im, threshold1));
t2 = uint8(imbinarize(original_im, threshold2));
% Sum the thresholded images, in this way we can select white parts in
% all the thresholded images
imf2 = preprocessed_im + t1 + t2;

% Dilate the filtered image to reduce black holes inside the wheel
imf2 = imbinarize(imf2);
se = strel('disk',2,8);
imf2 = imdilate(imf2,se);

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
figure(); imshow(original_im); title('Ellipsis of interest');
for k=1:length(ellipsis)
    drawEllipse(ellipsis(k),t);
end
clear s1; clear s2; clear t;

%% 2.1 Determine the ratio between the diameter of the wheel circles and the wheel-to-wheel distance 
% Get the conic matrices associated to the ellipsis found before. We can
% retrieve the equation of the ellipsis from their centers, axis lengths
% and orientation. Then, find the tangent lines by intersecting the dual
% conics

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
styles=['r-','m-','b-','y-']; widths=[1 2 2 1];
for k=1:length(tangent_lines)
    plotLine(tangent_lines(k,1),tangent_lines(k,2),1,[0 maxw], styles(k), widths(k));
end
clear styles; clear widths; clear k;

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
plot(upperwheelpoint1(1),upperwheelpoint1(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint1(1),lowerwheelpoint1(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(upperwheelpoint2(1),upperwheelpoint2(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint2(1),lowerwheelpoint2(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(middle_point(1),middle_point(2),'r^','MarkerSize',6,'HandleVisibility','off');
plot(vanish_pointz(1),vanish_pointz(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(vanish_pointy(1),vanish_pointy(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(1,1),ant_wheel_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(2,1),ant_wheel_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(1,1),post_wheel_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(2,1),post_wheel_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plotLine(line_vert_post_wheel(1),line_vert_post_wheel(2),line_vert_post_wheel(3),[0 maxw],'-',2);
plotLine(line_vert_ant_wheel(1),line_vert_ant_wheel(2),line_vert_ant_wheel(3),[0 maxw],'-',2);

% Get the image of the line at the infinity
inf_point1 = cross( lowertangentline, uppertangentline );
inf_point2 = cross( line_vert_post_wheel, line_vert_ant_wheel );
line_infinity = cross( inf_point1, inf_point2 );
plotLine(line_infinity(1),line_infinity(2),line_infinity(3),[0 maxw], '-', 2);
% Get the image of the circular points by intersecting the ellipsis and the
% line at the infinity
circ_point1 = intersectsconics(wheel_post_coeffs, [0 0 0 line_infinity]);
circ_point2 = intersectsconics(wheel_ant_coeffs, [0 0 0 line_infinity]);
circ_points = (circ_point1 + circ_point2) ./ 2;
legend('line1','line2','line3','line4', 'vertical1','vertical2', 'infinity');

hold off;

%% 2.2 Rectify the image using the image of the absolute conic.
% The image of the absolute conic can be computed by multiplying the images
% of the circular points:
% 
% $\omega_{\infty} = I*J' + J*I'$
% 
I = [circ_points(1,:)'; 1];
J = [circ_points(2,:)'; 1];
image_dual_conic = I*J.' + J*I.';
[~, DC, H] = svd(image_dual_conic);
normalization = sqrt(DC); normalization(3,3)=1;
%H_norm = inv(normalization*H);
H = normalization * H;
%H = inv(H);

% Find the inverse mapping of the four intersection points
a = H \ [post_wheel_points(2,:), 1]'; a = a/a(3);
b = H \ [post_wheel_points(1,:), 1]'; b = b/b(3);
c = H \ [ant_wheel_points(2,:), 1]'; c = c/c(3);
d = H \ [ant_wheel_points(1,:), 1]'; d = d/d(3);

% Compute the ratio distance - diameter from the rectified points
ratio = norm( a-c, 2) / norm( a-b, 2)

%figure();
%plot(a(1),a(2),'o'); plot(b(1),b(2),'x');
%plot(c(1),c(2),'h'); plot(d(1),d(2),'p');
%title('Inverse mapping of 4 wheel keypoints');

% 
%rectified_im = imwarp(original_im, projective2d(H'));

%figure(5);
% size_im = size(preprocessed_im);
% rectified_im = zeros(size_im);
% for i=1:size_im(1)
%     for j=1:size_im(2)
%         newpoint = H * [i; j; 1];
%         newpoint = newpoint / newpoint(3);
%         newpoint(1) = ceil(mod(newpoint(1), size_im(1))) +1;
%         newpoint(2) = ceil(mod(newpoint(2), size_im(2))) +1;
%         if newpoint(1) ~= 0 && newpoint(2) ~= 0
%             rectified_im(newpoint) = preprocessed_im(i,j);
%         end
%     end
% end
% for i=1:size_im(1)*4
%     for j=1:size_im(2)*4
%         index = [i/4; j/4; 1];
%         new_index = H\index;
%         new_index = ceil(new_index/new_index(3));
%         new_index(1) = mod(ceil(new_index(1)), size_im(1)) +1;
%         new_index(2) = mod(ceil(new_index(2)),size_im(2)) +1;
%         rectified_im(i,j) = preprocessed_im(new_index(1), new_index(2));
%     end
% end
%imshow(rectified_im);

%% Detect the features
% Get some features in the image using Harris method
figure();
[feat1_x,feat1_y] = harris(original_im, 2.7);
imshow(original_im), hold on, plot(feat1_y,feat1_x,'r+'); hold off;
% Take 4 feature points beloging to parallel lines in order to find another vanishing point
stop_point1 =  [320; 111]; %[306; 365];
stop_point2 =  [387; 117]; %[229; 341];
plate_point1 = [303; 393]; %[228; 365];
plate_point2 = [226; 367]; %[300; 389];

line_lights = cross( [stop_point1; 1], [stop_point2; 1] );
line_plate = cross( [plate_point1; 1], [plate_point2; 1] );

plotLine(line_lights(1), line_lights(2), line_lights(3), [-maxw maxw],'g-',2);
plotLine(line_plate(1), line_plate(2), line_plate(3), [-maxw maxw],'g-',2);
hold on;
vanish_pointx = cross( line_lights, line_plate );
vanish_pointx = vanish_pointx / vanish_pointx(3);
plot(vanish_pointx(1),vanish_pointx(2),'r*','MarkerSize',6,'HandleVisibility','off');

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
        %[circ_points(2,:) 1] * w * [circ_points(2,:) 1].';
      ];
res = solve(eqs);
fx = real(double(res.fx)); fx = fx(fx > 0); fx = fx(1);
fy = real(double(res.fy)); fy = fy(fy > 0); fy = fy(1);
u0 = real(double(res.u0)); u0 = u0(1);
v0 = real(double(res.v0)); v0 = v0(1);
clear K, clear eqs, clear res;
K = [fx, 0, u0; ...
     0, fy, v0; ...
     0, 0, 1]

% central_line = cross(cross(uppertangentline,lowertangentline),cross([tangent_lines(1,:),1],[tangent_lines(4,:),1]));
% pointspost = intersectsconics(wheel_post_coeffs, [0 0 0 central_line]);
% antspost = intersectsconics(wheel_ant_coeffs, [0 0 0 central_line]);

%% 3D Reconstruction
% Find the location of some image keypoints in the scene using
% backprojection. We fix the world reference frame at the center of the car
% plate, with the y-axis contained in the vertical symmetry plane.
% We can find the rotation matrix of the camera with respect to the world
% reference frame. We can exploit again the vanishing points, since
% (assuming that $v_x$ is an image point and $V_x$ is the direction of the
% x-axis in the scene):
%
% $$ v_x = P*V_x $$
% $$ v_x = K*\begin{pmatrix}R & t\end{pmatrix}*\begin{pmatrix}1 \\ 0 \\ 0 \\ 0\end{pmatrix} $$
% $$ v_x = K*\begin{pmatrix}r1 & r2 & r3 & t\end{pmatrix}*\begin{pmatrix}1 \\ 0 \\ 0 \\ 0\end{pmatrix} $$
% $$ v_x = K*r1 $$
% $$ r1 = K^{-1}*v_x $$
% 

% Compute the rows of R and normalize each one:
r1 = K \ vanish_pointx; r1 = r1/norm(r1);
r3 = K \ vanish_pointz; r3 = r3/norm(r3);
%r2 = K \ vanish_pointy; r2 = r2/norm(r2);
% Last column can be found by cross product by the other two
r2 = cross(r3,r1);
plate_center = [262; 365; 1];
t = K \ plate_center;
%t = [0; 0; 0];

R = [r1 r2 r3];
P = K*[R t];
M = P(:,1:3);

O = null(P); O = O/O(4);

% Get some symmetric features
plate_left = [218; 351; 1];
plate_right = [315; 383; 1];
% use stop points found before
stop_left = [stop_point1; 1];
stop_right = [stop_point2; 1];
light_right = [328; 262; 1];
light_left = [216; 236; 1];
%stem_center = [258; 259; 1];

% Find the backprojection rays of some symmetric points:
viewing_ray_platel = M \ plate_left;
viewing_ray_plater = M \ plate_right;
viewing_ray_lightl = M \ light_left;
viewing_ray_lightr = M \ light_right;
viewing_ray_stopl = M \ stop_left;
viewing_ray_stopr = M \ stop_right;

% Since we placed the origin of the world frame in the symmetry plane, its
% equation will be x=0. Knowning the viewing rays of the pairs of the projected
% image points that were symmetric in the scene, we can find the original
% 3D location by imposing that this points must be symmetric with respect
% to the symmetry plane x=0 (so, they must be opposite x coordinate sign).
% In addition, we also know that their y coordinates must be equal, again
% for their symmetric relation. This is equivalent to solve the system:
% 
% $$ \frac{x_1-x_O}{l_1} = \frac{y_1-y_O}{m_1} = \frac{z_1-z_O}{n_1} $$
%
% $$ \frac{x_2-x_O}{l_2} = \frac{y_2-y_O}{m_2} = \frac{z_2-z_O}{n_2} $$
%
% $$ x_1 = -x_2 $$
%
% $$ z_1 = z_2 $$
% 
% where $d_1 = (l_1, m_1, n_1)$ and $d_2 = (l_2, m_2, n_2)$ are the viewing rays
% passing though $P_1 = (x_1, y_1, z_1)$ and $P_2 = (x_2, y_2, z_2)$
% respectively, and $O = (x_O, y_O, z_O)$.
%
figure();
hold on;
plotPoint3([0,0,0],'wo',8);  % origin
% versors
plotLine3([0,0,0],[0.5,0,0],'r',2);
plotLine3([0,0,0],[0,0.5,0],'g',2);
plotLine3([0,0,0],[0,0,0.5],'b',2);
% viewpoint
plotPoint3(O,'m^',8);
% camera frame origin
plotPoint3(-t,'kh',1);

drawVRay(M,O,plate_left, 'c');
drawVRay(M,O,plate_right, 'c');

drawVRay(M,O,light_left, 'y');
drawVRay(M,O,light_right, 'y');

drawVRay(M,O,stop_left, 'r');
drawVRay(M,O,stop_right, 'r');

% planes z=0 and x=0
drawPlane([0; 0; 1; 0], 'k');
drawPlane([1; 0; 0; 0], 'c');

% Solve the system to find the pairs of symmetric features
[symm_point_platel,symm_point_plater] = findSymmetricFeatures(viewing_ray_platel, viewing_ray_plater, O);
[symm_point_lightl,symm_point_lightr] = findSymmetricFeatures(viewing_ray_lightl, viewing_ray_lightr, O);
[symm_point_stopl,symm_point_stopr] = findSymmetricFeatures(viewing_ray_stopl, viewing_ray_stopr, O);

% Plot the 3D locations of the symmetric features
plotPoint3(symm_point_platel,'mp',8);
plotPoint3(symm_point_plater,'mp',8);
plotPoint3(symm_point_lightl,'b^',6);
plotPoint3(symm_point_lightr,'b^',6);
plotPoint3(symm_point_stopl,'rs',6);
plotPoint3(symm_point_stopr,'rs',6);

hold off;
pbaspect([1,1,1]);
grid on;
title('3D Reconstruction')
xlabel('x'); ylabel('y'); zlabel('z');
legend('origin','x axis','y axis','z axis', 'viewpoint','cam frame', 'left plate ray','right plate ray','left light ray','right light ray','left stop ray','right stop ray', 'car back plane','symmetry plane', 'left plate','right plate','left light','right light','left stop','right stop');
set(gca, 'CameraPosition', [0,1,-1]);





