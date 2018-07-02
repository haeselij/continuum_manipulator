syms q_1 q_2 q_3 r_b v_1 v_2 v_3
% assume(q_1,'real');
% assume(q_2,'real');
% assume(q_3,'real');
% assume(v_1,'real');
% assume(v_2,'real');
% assume(v_3,'real');

theta = 2/3*(sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/r_b;
l_bar = (q_1 + q_2 + q_3)/3;
phi = atan((sqrt(3)*(q_3 - q_2))/(q_2 + q_3 - 2*q_1));
kappa = 2*(sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/(r_b*(q_1 + q_2 + q_3));

w_H_1h_1 = [cos(phi)^2*(cos(theta) - 1) + 1, sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)*sin(theta), (cos(phi)*(1 - cos(theta)))/kappa];
w_H_1h_2 = [sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)^2*(1 - cos(theta)) + cos(theta), sin(phi)*sin(theta), (sin(phi)*(1 - cos(theta)))/kappa];
w_H_1h_3 = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta), sin(theta)/kappa];
w_H_1h_4 = [0, 0, 0, 1];
w_H_1h = [w_H_1h_1; w_H_1h_2; w_H_1h_3; w_H_1h_4];

d_w_H_1h__dq1 =  diff(w_H_1h,q_1);
d_w_H_1h__dq2 =  diff(w_H_1h,q_2);
d_w_H_1h__dq3 =  diff(w_H_1h,q_3);

% d2w_H_1h__dq1dq1 = diff(d_w_H_1h__dq1,q_1)*[ 0; 0; 0; 1];
% d2w_H_1h__dq1dq2 = diff(d_w_H_1h__dq1,q_2)*[ 0; 0; 0; 1];
% d2w_H_1h__dq1dq3 = diff(d_w_H_1h__dq1,q_3)*[ 0; 0; 0; 1];
% d2w_H_1h__dq2dq2 = diff(d_w_H_1h__dq2,q_2)*[ 0; 0; 0; 1];
% d2w_H_1h__dq2dq3 = diff(d_w_H_1h__dq2,q_3)*[ 0; 0; 0; 1];
% d2w_H_1h__dq3dq3 = diff(d_w_H_1h__dq3,q_3)*[ 0; 0; 0; 1];

q_ =[q_1 q_2 q_3];
r_phi = q_1^2 + q_2^2 + q_3^2 - q_1*q_3 - q_2*q_3 - q_1*q_2;
r_k = -1/3*r_b*(q_1 + q_2 + q_3)^2*(r_phi)^(1/2);
dphi1 = 3^(1/2)*(q_2 - q_3)/(2*r_phi);
dphi2 = 3^(1/2)*(q_3 - q_1)/(2*r_phi);
dphi3 = 3^(1/2)*(q_1 - q_2)/(2*r_phi);
dkappa1 = (q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3)/r_k;
dkappa2 = (q_1^2 + q_3^2 - q_1*q_2 - q_2*q_3)/r_k;
d2phi1phi1 = 3^(1/2)/(2*r_phi^2)*(q_3 - q_2)*(2*q_1 - q_2 - q_3);
d2kappakappa1 = 3/(2*r_b*(q_1 + q_2 + q_3)^3*(r_phi^3)^(1/2))*((q_2 - q_3)^4 - 6*q_1*(q_2 + q_3)*(2/3*q_1^2 + q_2^2 + q_3^2) + 9*q_1^2*(q_2^2 + 2/3*q_2*q_3 + q_3^2) + 4*(q_2^4 + q_3^4));

for i=1:3
d_phi_dq(i) = -diff(phi,q_(i));
d_kappa_dq(i) = diff(kappa,q_(i));
d_l_bar_dq(i) = diff(l_bar,q_(i));
d2_phi_dq1q(i) = diff(phi,q_(i));
for k=1:3
d2_phi_ddq(i,k) = diff(d_phi_dq(i),q_(k));
d2_kappa_ddq(i,k) = diff(d_kappa_dq(i),q_(k));
d2_l_bar_ddq(i,k) = diff(d_l_bar_dq(i),q_(k));
end
end
% 
% d_k_dq1 = -diff(k_sys,q_1);
% d_k_dq2 = diff(k_sys,q_2);
% d_k_dq3 = diff(k_sys,q_3);
% 
% d_phi_dq1 = - diff(phi,q_1);
% d_phi_dq2 = - diff(phi,q_2);
% d_phi_dq3 = - diff(phi,q_3);
% d2_k_dq1dq1 =  diff(d_phi_dq1,q_1)
% d2_k_dq1dq2 = diff(d_k_dq1,q_2);
% d2_k_dq1dq3 = diff(d_k_dq1,q_3);
% d2_k_dq2dq2 = diff(d_k_dq2,q_2);
% d2_k_dq2dq3 = diff(d_k_dq2,q_3);
% d2_k_dq3dq3 = diff(d_k_dq3,q_3);

% d_kappa_dq1 = diff(kappa,q_1);
% d_kappa_dq2 = diff(kappa,q_2);
% d_l_bar_dq1 = diff(l_bar,q_1);
% subs(d_phi_dq3, [q_1, q_2,q_3],[0.30, 0.1, 0.3013])
% subs(dphi3, [q_1, q_2,q_3],[0.30, 0.10, 0.3013])
% r_b = 0.026;
% i =(subs(d2_kappa_ddq(1,1), [q_1, q_2,q_3],[0.1, 0.3014, 0.3013]))
%  k = ( subs(d2kappakappa1, [q_1, q_2,q_3],[0.1, 0.3014, 0.3013]))
% if i == k
%     n=10
% end
% roundn = @(x,n) round(x*10^n)./10^n;
% roundn(i/k,100)
