%% initialize value
syms  C_1 D_damp k_spline k
% 
% v_ = zeros(3, pressure_size); 

v_(:,1) = [0; 0; 0];

w_H_1h_ = zeros(4,4); 
dw_H_1h__dq1_ = zeros(4,4);
dw_H_1h__dq2_ = zeros(4,4);
dw_H_1h__dq3_ = zeros(4,4);

w_H_2dot_ = zeros(4,1);

% k_spline = zeros(length_size,1);
% C_1_ = zeros(length_size,1);
% D_damp_ = zeros(length_size,1);
% % 
% % coordinates_measured_cut(:,coordinate_size+1) = coordinates_measured_cut(:,coordinate_size);
% % coordinates_measured_cut(:,coordinate_size+2) = coordinates_measured_cut(:,coordinate_size);
% % 
% % length_time_converted = seconds(length_time_cut.');% transform_time_converted = seconds(transform_time_cut.');
% % r_tip_pos_time_table = timetable(length_time_converted, r_tip_pos_set_cut(1,:).', r_tip_pos_set_cut(2,:).', r_tip_pos_set_cut(3,:).');
% % r_tip_pos_time_table.Properties.VariableNames = {'x','y','z'};
% 
% length_time_converted = seconds(length_time);
% pressure_time_converted = seconds(pressure_time);
% 
% q_time_table = timetable(length_time_converted, q_(1,:).', q_(2,:).', q_(3,:).');
% pressure_time_table = timetable(pressure_time_converted, pressure_total(2,:).', pressure_total(7,:).', pressure_total(9,:).');
% 
% q_time_table.Properties.VariableNames = {'q1','q2','q3'};
% pressure_time_table.Properties.VariableNames = {'p1','p2','p3'};
% 
% sync_time_table = synchronize(q_time_table, pressure_time_table, 'regular', 'linear', 'SamplingRate', 2000);
% %%
%for i =  2000:pressure_size
   dt = 0.0005;
  
%   q_ = [sync_time_table.q1(i) ; sync_time_table.q2(i); sync_time_table.q3(i)];
%   q_next = [sync_time_table.q1(i+1) ; sync_time_table.q2(i+1); sync_time_table.q3(i+1)]

   %dt = (length_time(i)-length_time(i-1));
    v_(:,1) = [0;0;0];
    q_ = [0.295972500000003;0.310499999999995;0.297727500000002];
    q_next = q_;
    p_ =[0.125; 0; 0];
   [w_H_1h_, w_H_2dot_,dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, q_dot_, kappa_, l_bar_, phi_]= SetMatrices(dt, q_, q_next, v_(:,1), v_(:,1))
   %v_dot_ = [v_1_; v_2_; v_3_];
   [generalized_force_term_] = GetGeneralizedForceTerm(w_H_1h_, dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, q_(1), q_(2), q_(3), 10^5*p_(1), 10^5*p_(2), 10^5*p_(3), q_dot_);
   [bending_stiffness_derivative_] = GetBendingStiffnessDeriv(10^5*p_(1), 10^5*p_(2), 10^5*p_(3), kappa_*l_bar_, phi_,q_(1),q_(2),q_(3));
   [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, w_H_2dot_);
   [potential_energy_term_] = GetPotentialEnergyTerm(q_(1), q_(2), q_(3), dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, kappa_, l_bar_, phi_, bending_stiffness_derivative_);
   
   [ k_spline(1), C_1_(1) k] = LagrangeEquation(generalized_force_term_, kinetic_energy_term_, potential_energy_term_, q_dot_, v_(1,1), v_(2,1), v_(3,1), q_(1), q_(2), q_(3))
  
   
%end

5331766562847711/18014398509481984