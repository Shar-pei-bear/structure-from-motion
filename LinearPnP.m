function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
x = [x,ones(size(x,1),1)];
X = [X,ones(size(X,1),1)];
xc = K\x';
xc = xc';
u = xc(:,1);
v = xc(:,2);
w = xc(:,3);
A = [zeros(size(x,1),4), -w.*X, v.*X
     w.*X, zeros(size(x,1),4), -u.*X
     -v.*X, u.*X, zeros(size(x,1),4)];
[~,~,v] = svd(A);
P = reshape(v(:,end),4,3)';
R = P(:,1:3);
[u,d,v] = svd(R);
R = u*v';
t = P(:,4)./d(1,1);
if det(R) < 0
    R =- R;
    t =- t;
end
C = -R'*t;







