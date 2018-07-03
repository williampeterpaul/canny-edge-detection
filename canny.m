%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function emap = canny(img,sigma,fsize,tlow,thigh)
% Applies canny edge detection algorithm to the given image
%  img   - given image matrix
%  sigma - the value of sigma for the derivative of gaussian blur
%  fsize - given size of the gaussian blur kernal to be applied
%  lowt  - given lower threshold, used for hysteresis to avoid streaking
%  hight - given upper threshold, anything above which is an edge pixel
%
% Good parameters:
%     cameraman.tif :  canny('cameraman.tif',  1, 5, 10, 20)
%
% Note: Marker may need to convert outputs to uint8 to view other functions
% working
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function emap = canny(img,sigma,fsize,tlow,thigh)
% Input error checking, depending on amount of input arguments, sets
% default values, nargin = number of arguments input
if (nargin < 1)
    error('Please Enter Image Value');
elseif (nargin == 1)
    sigma = 1; fsize = 5; tlow = 10; thigh = 20;
elseif (nargin == 2)
    fsize = 5; tlow = 10; thigh = 20;
elseif (nargin == 3)
    tlow = 10; thigh = 20;
elseif (nargin == 4)
    thigh = 20;
end

% Checking for negative values
if(sigma < 0 || fsize < 0 || tlow < 0 || thigh < 0)
    error('Please Enter Positive Integer Values')
end

% If the threshold values are the wrong way round, correct them
if thigh < tlow
    temp = tlow;
    tlow = thigh;
    thigh = temp;
end

% Read in image
img = imread(img);

% Convert image to double (0-1)
img = double(img);

% Blur image using gaussian blur kernal function
blur = gaussBlurKernal(img, sigma, fsize);

% Get magnitude and direction using mag_dir function
[temp_mag,temp_dir] = mag_dir(blur);

% Use non-maximal suppression to get thin edges
edgeMap = thinEdges(temp_mag,temp_dir);

% Finally use hysteresis to avoid streaking in final image
emap = hystereis(edgeMap, tlow, thigh);

imshow(emap);

%% Make Gaussian Kernal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function creates and applies a gaussian kernal to an image to make blur
    function fimg = gaussBlurKernal(img, sigma, fsize)
        s = (fsize-1)/2;

        % Map out matrix to apply kernal
        [x,y] = meshgrid(-s:s, -s:s);
        k = exp(-(x.^2 + y.^2)/(2*sigma.^2));
        kernal = k/sum(sum(k));

        % Apply Gaussian Blur to image
        fimg = img;
        [x,y] = size(fimg);
        f = floor(fsize/2);

        % Using Padarray to pad the array correctly
        pimg = padarray(fimg,[f,f], 'replicate','both');

        % for-loop to iterate throughout size of image
        for i=1:x
            for j=1:y
                section = pimg(i:i+(fsize-1), j:j+(fsize-1));
                fimg(i,j) = kernal(:)' * section(:);
            end
        end
    end

%% Magnitude and Direction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applies appropriate equations to retieve matricies of magnitude and
% direction
    function [mag, dir] = mag_dir(fimg)
        fimg = padarray(fimg,[1,1], 'replicate','post');

        % Applies matrix values to get dy/dx
        dy = fimg(:,2:end) - fimg(:,1:end-1);
        dy = padarray(dy,[0,1], 'replicate','post');
        dx = fimg(2:end,:) - fimg(1:end-1,:);
        dx = padarray(dx,[1,0], 'replicate','post');

        dir = atan2d(dx,dy);
        mag = sqrt((dy.^2) + (dx.^2));

        % Applies roundUp function to the direction matrix
        dir = roundUp(dir);
    end

% Rounds up the angle to either 0, 45, 90, or 135
    function dir = roundUp(dir)
        dir(dir < 22.5 & dir >= -22.5) = 0;
        dir(dir < 67.5 & dir >= 22.5 | dir >= -67.5 & dir < -22.5) = 45;
        dir(dir < 112.5 & dir >= 67.5 | dir >= -112.5 & dir < -67.5) = 90;
        dir(dir < 157.5 & dir >= 112.5 | dir >= -157.5 & dir < -112.5) = 135;
        dir(dir >= 157.5 | dir < -112.5) = 0;
    end

%% Non-Maxmal Suppression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds definite edge pixels by comparing candidate pixel to its neighbours
    function edgeMap = thinEdges(mag, dir)
        [x, y] = size(mag);
        edgeMap = mag;
        % For-loop to iterate through size of image, with 1 pixel padding
        for i=2:x-1
            for j=2:y-1
                % Find co-ordinates of the next and prev pixels depening
                % on the direction of the gradient
                switch(dir(i,j))
                    case 0
                        prevPixel = mag(i,j-1);
                        nextPixel = mag(i,j+1);
                    case 45
                        prevPixel = mag(i-1,j-1);
                        nextPixel = mag(i+1,j+1);
                    case 90
                        prevPixel = mag(i-1,j);
                        nextPixel = mag(i+1,j);
                    case 135
                        prevPixel = mag(i+1,j-1);
                        nextPixel = mag(i-1,j+1);
                    otherwise
                        error('Invalid Direction');
                end
                % If either neighbours are bigger, turn off candidate pixel
                if(mag(i,j) <= prevPixel || mag(i,j) < nextPixel)
                    edgeMap(i,j) = 0;
                end
            end
        end
    end

%% Hysterisis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compares candidate pixel to thresholds to ensure it is an edge pixel, if
% it is not above thigh, but it is above tlow, see if its neighbours are
% definite edge pixels, if they are, candidate becomes definite
    function emap = hystereis(edgeMap, tlow, thigh)

        % Get all pixels in edge map that are above the lower threshold
        aboveLower = edgeMap > tlow;
        % Get Co-ordinates of the pixels that are above the upper threshold
        [aboveUpper_X, aboveUpper_Y] = find(edgeMap > thigh);
        % Return emap of all points in aboveLower which are connected to a
        % point above the upper threshold
        emap = bwselect(aboveLower, aboveUpper_Y, aboveUpper_X);
    end

end
