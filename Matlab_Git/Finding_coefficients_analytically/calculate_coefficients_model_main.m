%% initialize value
syms  C_1 D_damp 


v_ = zeros(3, pressure_size); 
v_(:,1) = [0; 0; 0];

w_H_1h_ = zeros(4,4); 
dw_H_1h__dq1_ = zeros(4,4);
dw_H_1h__dq2_ = zeros(4,4);
dw_H_1h__dq3_ = zeros(4,4);

w_H_2dot_ = zeros(4,1);

k_ = zeros(length_size,1);
C_1_ = zeros(length_size,1);
D_damp_ = zeros(length_size,1);
% 
% coordinates_measured_cut(:,coordinate_size+1) = coordinates_measured_cut(:,coordinate_size);
% coordinates_measured_cut(:,coordinate_size+2) = coordinates_measured_cut(:,coordinate_size);
% 
length_time_converted = seconds(length_time_cut.');% transform_time_converted = seconds(transform_time_cut.');
r_tip_pos_time_table = timetable(length_time_converted, r_tip_pos_set_cut(1,:).', r_tip_pos_set_cut(2,:).', r_tip_pos_set_cut(3,:).');
r_tip_pos_time_table.Properties.VariableNames = {'x','y','z'};

length_time_converted = seconds(length_time.');
pressure_time_converted = seconds(pressure_time.');

q_time_table = timetable(length_time, q_(1,:).', q_(2,:).', q_(3,:).');
pressure_time_table = timetable(pressure_time, pressure_total(2,:).', pressure_total(7,:).', pressure_total(9,:).');

r_tip_pos_time_table.Properties.VariableNames = {'x','y','z'};

coordinates_measured_time_table = timetable(transform_time_converted.', coordinates_measured_cut(1,:).', coordinates_measured_cut(2,:).', coordinates_measured_cut(3,:).'); 
coordinates_measured_time_table.Properties.VariableNames = {'x','y','z'};

%sync_time_table = synchronize(coordinates_measured_time_table, r_tip_pos_time_table, 'union', 'linear');
sync_time_table = synchronize(q_measured_time_table, pressure_time_table, 'regular', 'linear', 'SamplingRate', 5000);

for i =  2000:pressure_size
    
   v_(:,i) = (q_(:,i)-q_(:,i-1))/dt;
   dt = (length_time(i)-length_time(i-1));
   
   [w_H_1h_, w_H_2dot_,dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, q_dot_, kappa_, l_bar_, phi_]= SetMatrices(dt, q_(:,i), q_(:,i+1), v_(:,i), v_(:,i+1))
   %v_dot_ = [v_1_; v_2_; v_3_];
   [generalized_force_term_] = GetGeneralizedForceTerm(w_H_1h_, dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, q_(1,i), q_(2,i), q_(3,i), 10^5*pressure_total(2,i), 10^5*pressure_total(7,i), 10^5*pressure_total(9,i), q_dot_);
   [bending_stiffness_derivative_] = GetBendingStiffnessDeriv(10^5*pressure_total(2,i), 10^5*pressure_total(7,i), 10^5*pressure_total(9,i), kappa_*l_bar_, phi_,q_(1,i),q_(2,i),q_(3,i));
   [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, w_H_2dot_);
   [potential_energy_term_] = GetPotentialEnergyTerm(q_(1,i), q_(2,i), q_(3,i), dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, kappa_, l_bar_, phi_, bending_stiffness_derivative_);
   
   [ C_1_(i-1999,1), D_damp_(i-1999,1)] = LagrangeEquation(generalized_force_term_, kinetic_energy_term_, potential_energy_term_, q_dot_, v_(1,i), v_(2,i), v_(3,i), q_(1,i), q_(2,i), q_(3,i));
  
   
end