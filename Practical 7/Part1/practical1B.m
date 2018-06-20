function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');
figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 

hom1 = calcBestHomography(pts1, pts2);
hom2 = calcBestHomography(pts1b, pts3);

oneX = size(im1,1);
oneY = size(im1,2);
twoX = size(im2,1);
twoY = size(im2,2);
threeX = size(im3,1);
threeY = size(im3,2);

% Picture 2 to 1

% Transform pixel position and normalize
transform = hom1*[reshape(meshgrid(1:oneY,1:oneX),[1,oneX*oneY]);repmat(1:oneX,[1,oneY]);ones(1,oneX*oneY)];
transform(1:2,:) = ceil(transform(1:2,:)./repmat(transform(3,:),[2,1]));

% Transform position to new image
for y=1:oneY
    for x=1:oneX
        updatedX = transform(2,x+(y-1)*oneX);
        updatedY = transform(1,x+(y-1)*oneX);
        if(updatedX > 0 && updatedX < twoX ...
        && updatedY > 0 && updatedY < twoY)
            im1(x,y,:) = im2(updatedX,updatedY,:);
        end
    end
end


% Picture 3 to 1

% Transform pixel position and normalize
transform = hom2*[reshape(meshgrid(1:oneY,1:oneX),[1,oneX*oneY]);repmat(1:oneX,[1,oneY]);ones(1,oneX*oneY)];
transform(1:2,:) = ceil(transform(1:2,:)./repmat(transform(3,:),[2,1]));

% Transform position to new image
for y=1:oneY
    for x=1:oneX
        updatedX = transform(2,x+(y-1)*oneX);
        updatedY = transform(1,x+(y-1)*oneX);
        if(updatedX > 0 && updatedX < threeX ...
        && updatedY > 0 && updatedY < threeY)
            im1(x,y,:) = im3(updatedX,updatedY,:);
        end
    end
end

figure; image(uint8(im1)); axis off;

% Function same as in practical1.m
function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding matchin in 
%pts2Cart

%**** TO DO ****;

% turn points to homogeneous
pts1Cart = [pts1Cart; ones(1,5)];
pts2Cart = [pts2Cart; ones(1,5)];

% construct 10 x 9 matrix
A = zeros(10,9);
for count = 1:size(pts1Cart,2)
    A(2*count,:) = [pts1Cart(1,count),pts1Cart(2,count),1,0,0,0,-pts2Cart(1,count)*pts1Cart(1,count),-pts2Cart(1,count)*pts1Cart(2,count),-pts2Cart(1,count)];
    A((2*count)-1,:) = [0,0,0,-pts1Cart(1,count),-pts1Cart(2,count),-1,pts2Cart(2,count)*pts1Cart(1,count),pts2Cart(2,count)*pts1Cart(2,count),pts2Cart(2,count)];
end

%solve Ah = 0
h = solveAXEqualsZero(A);

%reshape h into the matrix H
H = (reshape(h,[3,3]))';


% Function same as in from practical1.m
function x = solveAXEqualsZero(A);

[U,S,V] = svd(A);
x = V(:,9);