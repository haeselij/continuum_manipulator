%% initialize value
syms q_1_n q_2_n q_3_n v_1_n v_2_n v_3_n 

% random value for debuging TODO
% pressure_size = 100;
% pressure_total= zeros(12,pressure_size);
% pressure_total(2,:) = 0;
% pressure_total(7,:) = 0;
% pressure_total(9,:) = 0;
% t = 0;

q_container = zeros(3, pressure_size); 
v_container = zeros(3, pressure_size); 
w_H_1h_ = zeros(4,4); % coordinatetransform matrix head to worldframe

dw_H_1h__dq1_ = zeros(4,4); %first derrivative of H_1h after q_1
dw_H_1h__dq2_ = zeros(4,4);
dw_H_1h__dq3_ = zeros(4,4);

w_H_2dot_ = zeros(4,1);

v_1_ = 0.000;
v_2_ = 0.000;
v_3_ = 0.000;

q_1_ = 0.3014; 
q_2_ = 0.3014;
q_3_ = 0.3014;

generalized_force_term_container = zeros(3, pressure_size);
bending_stiffness_derivative_container = zeros(3, pressure_size);
bending_stiffness_derivative_container = zeros(3, pressure_size);
kinetic_energy_term_container = zeros(3, pressure_size);
potential_energy_term_container = zeros(3, pressure_size);

for i =  2:pressure_size
   % random dt TODO
   %dt = 0.001;
   dt = (pressure_time(i)-pressure_time(i-1));
   
   q_container(1,i-1) = q_1_;
   q_container(2,i-1) = q_2_; 
   q_container(3,i-1) = q_3_;  
   
   v_container(1,i-1) = v_1_;
   v_container(2,i-1) = v_2_;
   v_container(3,i-1) = v_3_;
  
   [w_H_1h_, w_H_2dot_,dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, q_dot_, kappa_, l_bar_, phi_]= SetMatrices(dt, q_1_, q_2_, q_3_, v_1_, v_2_, v_3_);
   v_dot_ = [v_1_; v_2_; v_3_];
   [generalized_force_term_] = GetGeneralizedForceTerm(w_H_1h_, dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, q_1_, q_2_, q_3_, 10^5*pressure_total(2,i-1), 10^5*pressure_total(7,i-1), 10^5*pressure_total(9,i-1), v_dot_)
   [bending_stiffness_derivative_] = GetBendingStiffnessDeriv(10^5*pressure_total(2,i-1), 10^5*pressure_total(7,i-1), 10^5*pressure_total(9,i-1), kappa_*l_bar_, phi_,q_1_,q_2_,q_3_);
   [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_, w_H_2dot_);
   [potential_energy_term_] = GetPotentialEnergyTerm(real(q_1_), real(q_2_), real(q_3_), dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, kappa_, l_bar_, phi_, bending_stiffness_derivative_)
   
   [q_1_, q_2_, q_3_, v_1_, v_2_, v_3_] = LagrangeEquation(generalized_force_term_, kinetic_energy_term_, potential_energy_term_, q_dot_, v_1_, v_2_, v_3_, q_1_, q_2_, q_3_)
  
   generalized_force_term_container(:,i-1) = generalized_force_term_;
   bending_stiffness_derivative_container(:,i) = bending_stiffness_derivative_;
   kinetic_energy_term_num = GetNumKineticEnergy(q_1_, q_2_, q_3_, v_1_, v_2_, v_3_,kinetic_energy_term_);
   kinetic_energy_term_container(:,i) = kinetic_energy_term_num;
   potential_energy_term_container(:,i) = potential_energy_term_;
   
   %p = t*0.15
   %pressure_total(2,i) = t*0.15; %let pressure increase linearly
   
   %t = t + dt
end