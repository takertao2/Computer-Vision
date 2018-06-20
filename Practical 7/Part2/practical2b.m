function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.

TEst = estimatePlanePose(xImCart,XCart,K);


%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];

%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points

projectionMult = [XWireFrameCart;ones(1,size(XWireFrameCart,2))]
drawData = [K, zeros(3,1)]*TEst*projectionMult;
drawQuot = repmat(drawData(3,:),3,1)
draw = drawData(1:3,:)./drawQuot;

plot([draw(1,1), draw(1,2)],[draw(2,1), draw(2,2)])
plot([draw(1,1), draw(1,4)],[draw(2,1), draw(2,4)])
plot([draw(1,1), draw(1,5)],[draw(2,1), draw(2,5)])
plot([draw(1,2), draw(1,3)],[draw(2,2), draw(2,3)])
plot([draw(1,2), draw(1,6)],[draw(2,2), draw(2,6)])
plot([draw(1,3), draw(1,4)],[draw(2,3), draw(2,4)])
plot([draw(1,3), draw(1,7)],[draw(2,3), draw(2,7)])
plot([draw(1,4), draw(1,8)],[draw(2,4), draw(2,8)])
plot([draw(1,5), draw(1,6)],[draw(2,5), draw(2,6)])
plot([draw(1,6), draw(1,7)],[draw(2,6), draw(2,7)])
plot([draw(1,5), draw(1,8)],[draw(2,5), draw(2,8)])
plot([draw(1,7), draw(1,8)],[draw(2,7), draw(2,8)])
plot(draw(1,:),draw(2,:),'x')


%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?

%Function same as in practical2.m
function T = estimatePlanePose(xImCart,XCart,K)

%TO DO Convert Cartesian image points xImCart to homogeneous representation
%xImHom

xImHom = [xImCart;ones(1,size(xImCart,2))];

%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom

xCamHom = inv(K)*xImHom; %pre-multiply

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.

H = calcBestHomography(XCart,xCamHom);

%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD

[U,S,V] = svd(H(:,1:2));
R = U*[1 0; 0 1; 0 0]*V';

%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns

crossProd = cross(R(:,1),R(:,2));
R = [R, crossProd];

%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
if(det(R) <= 0)
    R(:,3) = -1*R(:,3);
end

%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H

scalingSum = sum(H(:,1:2)./R(:,1:2))
scalingQuot = sum(scalingSum)/6;
t = H(:,3)/scalingQuot;

%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.

if(t(3) <= 0)
    t = -1*t;
    R(:,1:2) = -1*R(:,1:2);
end

%assemble transformation into matrix form
%T  = [R t;0 0 0 1];

T  = [R t;0 0 0 1];


% Function same as in practical1B.m
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
    % Different equation for odd and even, repeats (see slides)
    A(2*count,:) = [pts1Cart(1,count),pts1Cart(2,count),1,0,0,0,-pts2Cart(1,count)*pts1Cart(1,count),-pts2Cart(1,count)*pts1Cart(2,count),-pts2Cart(1,count)];
    A((2*count)-1,:) = [0,0,0,-pts1Cart(1,count),-pts1Cart(2,count),-1,pts2Cart(2,count)*pts1Cart(1,count),pts2Cart(2,count)*pts1Cart(2,count),pts2Cart(2,count)];
end

%solve Ah = 0
h = solveAXEqualsZero(A);

%reshape h into the matrix H
H = (reshape(h,[3,3]))';


% Function same as in from practical1B.m
function x = solveAXEqualsZero(A);

[U,S,V] = svd(A);
x = V(:,9);

