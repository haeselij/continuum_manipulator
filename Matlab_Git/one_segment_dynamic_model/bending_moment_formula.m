%% Derivation of Bending Moment
syms D_1 D_2  C_01 C_02 C_1 C_2 q__1 q__2 q__3 r_b kappa theta_1 theta_2 theta
assume(q__1,'real');
assume(q__2,'real');
assume(q__3,'real');

theta = 2/3*(sqrt(q__1^2 + q__2^2 + q__3^2 - q__1*q__2 - q__1*q__3 - q__2*q__3))/r_b;
l_bar = (q__1 + q__2 + q__3)/3;
phi = atan2((sqrt(3)*(q__3 - q__2)),(q__2 + q__3 - 2*q__1));
k_sys = [phi; kappa; l_bar];

%theta_2 = theta - theta_1; 
%sine rule
eqn_1 = theta_1/sin(pi/3- phi) == theta_2/sin(phi);
eqn_2 = theta_1/sin(pi/3- phi) == theta/sin(2/3*pi);
[theta_1, theta_2] = solve([eqn_1 eqn_2], [theta_1 theta_2])

%cosine rule
D_2 = 1 - D_1; 

M_1 = theta*C_1;

M_2 = theta_2*C_2;

M_bend = M_1 + M_2;

d_M_bend_dq__1 = diff(M_1, q__1);
d_M_bend_dq__2 = diff(M_2, q__2);
d_M_bend_dq__3 = diff(M_2, q__3);

% q__1 = 0.3014;
% q__2 = 0.3014;
% q__3 = 0.3014;

M_1 = d_M_bend_dq__1



% M_2 = theta*cos(3/(2*pi)*phi)*C_2;
% 
% 
% M_bend = M_1 + M_2;
% 
% % for i=1:3
% d_M_bend_dq1 = diff(M_bend, q__1);
% d_M_bend_dq2 = diff(M_bend, q__2)
% d_M_bend_dq3 = diff(M_bend, q__1);

% end

%% getting the fucking derivative of M_Bend
% huge mess of cases
%case 1 

% if phi < pi/3 && phi >= -pi/3
%     d_M_bend_dq1 = 0;
%     if phi < pi/3 && phi >= 0
%         q_1 = q__2;
%         q_2 = q__3;
%         q_3 = q__1;
% 
%     else
%         q_1 = q__2;
%         q_2 = q__1;
%         q_3 = q__3;
%     end         
%         d_M_bend_dq2 = diff(M_bend, q_2);
%         d_M_bend_dq3 = diff(M_bend, q_3);
%         
% else if phi < 2*pi/3 &&  phi >= pi/3
%         d_M_bend_dq2 = 0;
%     end 
