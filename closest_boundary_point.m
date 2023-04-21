% computes the convex hull of a group of blobs, and specifies a point outside that convex hull.
% Then finds the pixel on the convex hull closest to that external point,
% computes the distance between them and draws a line between them.
clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.
close all;
clear all;
format long g;
format compact;
fontSize = 20;

grayImage = peaks(300);
% Display the original gray scale image.
subplot(2, 2, 1);
imshow(grayImage, []);
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
% Give a name to the title bar.
set(gcf,'name','Demo by ImageAnalyst','numbertitle','off') 

binaryImage = grayImage > 0.95;
% Display the image.
subplot(2, 2, 2);
imshow(binaryImage, []);
title('Binary Image', 'FontSize', fontSize);

% Compute the convex hull of the group of 2 blobs:
chImage = bwconvhull(binaryImage);
% Display the image.
subplot(2, 2, 3);
imshow(chImage, []);
axis on;
title('Convex Hull Image', 'FontSize', fontSize);

% Mark a point outside the convex hull
yRow = 50;
xCol = 150;
hold on;
plot(xCol, yRow,'r+', 'MarkerSize', 15, 'LineWidth', 2);

% Calculate the distance transform
[edtImage, indexImage] = bwdist(chImage);
% Display the image.
subplot(2, 2, 4);
imshow(edtImage, []);
axis on;
title('Euclidean Distance Transform', 'FontSize', fontSize);

% Find the distance from (150,50) to the point nearest to it on the convex hull
theDistance = edtImage(yRow, xCol)

% Plot the line from the point to the point on the boundary
[row, col] = ind2sub(size(grayImage), indexImage(yRow, xCol));
hold on;
line([xCol, col], [yRow, row], 'Color', 'r', 'LineWidth', 2);
plot(xCol, yRow,'r+', 'MarkerSize', 15, 'LineWidth', 2);
% Get the boundary
boundaries = bwboundaries(chImage);
xb = boundaries{1}(:, 2);
yb = boundaries{1}(:, 1);
plot(xb, yb, 'g-');

% Tell user the results.
message = sprintf('The distance = %.3f pixels.', theDistance);
uiwait(helpdlg(message));

stats = regionprops(chImage,"MajorAxisLength", "MinorAxisLength", "EquivDiameter");