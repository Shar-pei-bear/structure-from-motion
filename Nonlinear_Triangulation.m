function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 
syms x y z 
x_3d = [x;y;z;1];

P1 = K * [R1 -R1*C1];
P2 = K * [R2 -R2*C2];
P3 = K * [R3 -R3*C3];

x_c1_h = P1*x_3d;
x_c1   = x_c1_h(1:2)./x_c1_h(3);
x_c2_h = P2*x_3d;
x_c2   = x_c2_h(1:2)./x_c2_h(3);
x_c3_h = P3*x_3d;
x_c3   = x_c3_h(1:2)./x_c3_h(3);

x_c = [x_c1;x_c2;x_c3];
J = jacobian(x_c,[x,y,z]);
X = X0;
b = [x1,x2,x3];
for i = 1:size(x1,1)
    step_length = 1;
    iteration = 0;
    x_c_r = eval(subs(x_c,[x,y,z],X(i,:)));
    while step_length > 1e-05 && iteration < 100
        J_r = eval(subs(J,[x,y,z],X(i,:)));  
        delta_x = (J_r'*J_r)\J_r'*(b(i,:)'-x_c_r);
        X(i,:) = X(i,:) + delta_x';
        x_c_r_old = x_c_r;
        x_c_r     = eval(subs(x_c,[x,y,z],X(i,:)));
        step_length = norm(x_c_r - x_c_r_old)/norm(x_c_r_old);
        iteration = iteration + 1;
    end
end