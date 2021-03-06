% In this script the second derivative of the transform matrix is derived.
% The result is then copied into the SetMatrices function directly.

syms  phi kappa l_bar

theta = kappa*l_bar;
q = [phi; kappa; l_bar];
w_H_1h_1 = [cos(phi)^2*(cos(theta) - 1) + 1, sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)*sin(theta), (cos(phi)*(1 - cos(theta)))/kappa];
w_H_1h_2 = [sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)^2*(1 - cos(theta)) + cos(theta), sin(phi)*sin(theta), (sin(phi)*(1 - cos(theta)))/kappa];
w_H_1h_3 = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta), sin(theta)/kappa - 0.3014];
w_H_1h_4 = [0, 0, 0, 1];
w_H_1h = [w_H_1h_1; w_H_1h_2; w_H_1h_3; w_H_1h_4];

d_w_H_1h_dkappa = diff(w_H_1h,kappa);
d_w_H_1h_dphi = diff(w_H_1h,phi);
d_w_H_1h_dl_bar = diff(w_H_1h,l_bar);

d_w_H_1h__dq1 =  diff(w_H_1h,q_1);
d_w_H_1h__dq2 =  diff(w_H_1h,q_2);
d_w_H_1h__dq3 =  diff(w_H_1h,q_3);

d2_w_H_1h_dkappadkappa = diff(d_w_H_1h_dkappa, kappa);
d2_w_H_1h_dphidphi = diff(d_w_H_1h_dphi , phi);
d2_w_h_1h_dl_bard_lbar = diff(d_w_H_1h_dl_bar, l_bar);

d2_w_H_1h_dphidkappa = diff(d_w_H_1h_dphi, kappa);
d2_w_H_1h_dphidl_bar = diff(d_w_H_1h_dphi, l_bar);
d2_w_H_1h_dkappadl_bar = diff(d_w_H_1h_dkappa, l_bar);

for i=1:3
    for k=1:3
        
        d_2_H_1h_dkdk_dd(i,k,:,:) = d2_w_H_1h_dphidphi*d_phi_dq(i)*d_phi_dq(k) + d2_w_H_1h_dphidkappa*d_phi_dq(i)*d_kappa_dq(k) + d2_w_H_1h_dphidl_bar*d_phi_dq(i)*d_l_bar_dq(k) + d2_w_H_1h_dphidkappa*d_kappa_dq(i)*d_phi_dq(k) + d2_w_H_1h_dkappadkappa*d_kappa_dq(i)*d_kappa_dq(k) + d2_w_H_1h_dkappadl_bar*d_kappa_dq(i)*d_l_bar_dq(k)   + d2_w_H_1h_dphidl_bar*d_l_bar_dq(i)*d_phi_dq(k) + d2_w_H_1h_dkappadl_bar*d_l_bar_dq(i)*d_kappa_dq(k) + d2_w_h_1h_dl_bard_lbar*d_l_bar_dq(i)*d_l_bar_dq(k);
        
        dH_1h_dk_dd(i,k,:,:) = d_w_H_1h_dphi*d2_phi_ddq(i,k) + d_w_H_1h_dkappa*d2_kappa_ddq(i,k) + d_w_H_1h_dl_bar*d2_l_bar_ddq(i,k);
        
        d2_H_1h_dq(i,k,:,:) = d_2_H_1h_dkdk_dd(i,k,:,:) + dH_1h_dk_dd(i,k,:,:);
        
    end
end

for i=1:4   
    
    d2w_H_1h__dq1dq1(i,1) = d2_H_1h_dq(1,1,i,4);
    d2w_H_1h__dq1dq2(i,1) = d2_H_1h_dq(1,2,i,4);
    d2w_H_1h__dq1dq3(i,1) = d2_H_1h_dq(1,3,i,4);
    d2w_H_1h__dq2dq2(i,1) = d2_H_1h_dq(2,2,i,4);
    d2w_H_1h__dq2dq3(i,1) = d2_H_1h_dq(2,3,i,4);
    d2w_H_1h__dq3dq3(i,1) = d2_H_1h_dq(3,3,i,4);
    
end
