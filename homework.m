%%  IMAGE ANALYSIS & COMPUTER VISION HOMEWORK
%   Federico Parroni 10636497
%   

close all
clear all
clc

%% 1. Ellipsis detection
% Load the original image, resize to 1/4 and convert it
% to grayscale
original_im = rgb2gray(imresize(imread('Image1.jpg'),0.25));
maxx = length(original_im) * 1.2;

figure(), imshow(original_im,[]);
title('Original image');

%% 1.1 Preprocess the original image
% Enhance the contrast of the image by remapping pixels intensities
curves = uint8([zeros(1,60) linspace(0,255,16) repelem(255,180)]);
figure();
plot(0:255, curves);
title('Preprocessing (curves)');
xlabel('original pixel intensity'); ylabel('transformed pixel intensity');

preprocessed_im = intlut(original_im, curves);
figure();
imshow(preprocessed_im);
title('Preprocessed image');

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

figure(), imshow(imf1);
figure(), imshow(imf2);


%% 1.2 Detect the ellipsis in the preprocessed images
% Get ellipsis using regionprops function
s1 = regionprops(imf1, {'Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity'});
s2 = regionprops(imf2, {'Centroid','MajorAxisLength','MinorAxisLength','Orientation','Eccentricity'});

% Filter all the ellipses detected before to exclude bad candidates
s1 = s1(arrayfun(@(x) x.MajorAxisLength < 120 && x.MajorAxisLength > 20 ...
                   && x.MinorAxisLength < 80 && x.MajorAxisLength > 20 ...
                   && (x.Orientation < -80 && x.Orientation > -90 ...
                   || x.Orientation < 78 && x.Orientation > 72), s1));

s2 = s2(arrayfun(@(x) x.MajorAxisLength < 120 && x.MajorAxisLength > 110, s2));

% Draw the features
t = linspace(0,2*pi,60);
figure(), imshow(imf1);
for k=1:length(s1)
    drawEllipse(s1(k),t);
end
figure(), imshow(imf2);
for k=1:length(s2)
    drawEllipse(s2(k),t);
end

% Combine both results
ellipsis = [s1; s2];
figure();
title('Ellipsis'), imshow(original_im);
for k=1:length(ellipsis)
    drawEllipse(ellipsis(k),t);
end

% Get the conic matrices of the two wheels
wheel_ant = ellipsis(2);
wheel_post = ellipsis(3);
wheel_ant_coeffs = ellipseCoeff(wheel_ant);
wheel_post_coeffs = ellipseCoeff(wheel_post);
dual_wheel_ant_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(wheel_ant_coeffs)));
dual_wheel_post_coeffs = conicCoeffFromMatrix(inv(conicMatrixFromCoeff(wheel_post_coeffs)));

tangent_lines = intersectsconics(dual_wheel_ant_coeffs,dual_wheel_post_coeffs);
% Plot the 4 tangent lines
styles=['r-','m-','b-','y-'];
widths=[1 2 2 1];
for k=1:length(tangent_lines)
    plotLine(tangent_lines(k,1),tangent_lines(k,2),1,[0 maxx], styles(k), widths(k));
end
clear styles, clear widths, clear k;

%% 2. Determine the ratio between the diameter of the wheel circles and the wheel-to-wheel distance 
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
vanish_point1 = cross(uppertangentline, lowertangentline);
vanish_point1 = vanish_point1 / vanish_point1(3);
vanish_point2 = cross(line_vert_post_wheel, line_vert_ant_wheel);
vanish_point2 = vanish_point2 / vanish_point2(3);
% Get the 4 intersection points
middle_point = cross( [tangent_lines(1,:), 1], [tangent_lines(4,:), 1] );
middle_point = middle_point / middle_point(3);
middle_line = cross( vanish_point1, middle_point );
post_wheel_points = intersectsconics(wheel_post_coeffs, [0 0 0 middle_line]);
ant_wheel_points = intersectsconics(wheel_ant_coeffs, [0 0 0 middle_line]);

% Plot all tangency and vanishing points and lines
hold on;
plot(upperwheelpoint1(1),upperwheelpoint1(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint1(1),lowerwheelpoint1(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(upperwheelpoint2(1),upperwheelpoint2(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(lowerwheelpoint2(1),lowerwheelpoint2(2),'m^','MarkerSize',6,'HandleVisibility','off');
plot(middle_point(1),middle_point(2),'r^','MarkerSize',6,'HandleVisibility','off');
plot(vanish_point1(1),vanish_point1(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(vanish_point2(1),vanish_point2(2),'b*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(1,1),ant_wheel_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(ant_wheel_points(2,1),ant_wheel_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(1,1),post_wheel_points(1,2),'r*','MarkerSize',6,'HandleVisibility','off');
plot(post_wheel_points(2,1),post_wheel_points(2,2),'r*','MarkerSize',6,'HandleVisibility','off');
plotLine(line_vert_post_wheel(1),line_vert_post_wheel(2),line_vert_post_wheel(3),[0 maxx],'-',2);
plotLine(line_vert_ant_wheel(1),line_vert_ant_wheel(2),line_vert_ant_wheel(3),[0 maxx],'-',2);

% Get the image of the line at the infinity
inf_point1 = cross( lowertangentline, uppertangentline );
inf_point2 = cross( line_vert_post_wheel, line_vert_ant_wheel );
line_infinity = cross( inf_point1, inf_point2 );
plotLine(line_infinity(1),line_infinity(2),line_infinity(3),[0 maxx], '-', 2);
% Get the image of the circular points by intersecting the ellipsis and the
% line at the infinity
circ_point1 = intersectsconics(wheel_post_coeffs, [0 0 0 line_infinity]);
circ_point2 = intersectsconics(wheel_ant_coeffs, [0 0 0 line_infinity]);
circ_points = (circ_point1 + circ_point2) ./ 2;
legend('line1','line2','line3','line4', 'vertical1','vertical2', 'infinity');

hold off;

%% Rectify the image using the images of the circular points
I = [circ_points(1,:)'; 1];
J = [circ_points(2,:)'; 1];
image_dual_conic = I*J.' + J*I.';
[~, DC, H] = svd(image_dual_conic);
%normalization = sqrt(DC);
%H = inv(normalization*H);
%H = normalization * H;
H = inv(H);

% Find the inverse mapping of the four intersection points
a = H * [post_wheel_points(2,:), 1]';
b = H * [post_wheel_points(1,:), 1]';
c = H * [ant_wheel_points(2,:), 1]';
d = H * [ant_wheel_points(1,:), 1]';

% Compute the ratio distance - diameter from the rectified points
ratio = sqrt(sum( (a/a(3)-c/c(3)).^2 )) / sqrt(sum( (a/a(3)-b/b(3)).^2 ));
% 
% rectified_im = imwarp(original_im, projective2d(H));

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
% Get features using Harris method
%figure(5);
% try also with original_im, preprocessed_im
%[feat1_x,feat1_y] = harris(imf1, 3);
%imshow(imf1), hold on, plot(feat1_y,feat1_x,'r+'), hold off;


% central_line = cross(cross(uppertangentline,lowertangentline),cross([tangent_lines(1,:),1],[tangent_lines(4,:),1]));
% pointspost = intersectsconics(wheel_post_coeffs, [0 0 0 central_line]);
% antspost = intersectsconics(wheel_ant_coeffs, [0 0 0 central_line]);



