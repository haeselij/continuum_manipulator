%% initialize values 
q_container= zeros(3, pressure_size);  


w_H_1h_ = zeros(4,4); % coordinatetransform matrix head to worldframe

dw_H_1h__dq1_ = zeros(4,4); %first derrivative of H_1h after seg.q_1
dw_H_1h__dq2_ = zeros(4,4);
dw_H_1h__dq3_ = zeros(4,4);

w_H_2dot_ = zeros(4,1);

v_1_ = 0;
v_2_ = 0;
v_3_ = 0;

syms q_1_ q_2_ q_3_;  

for i = 2:10
   dt = pressure_time(i)-pressure_time(i-1);
  
   [w_H_1h_, w_H_2dot_,dw_H_1h__dq1_, dw_H_1h__dq2_, dw_H_1h__dq3_]= SetMatrices(dt, q_1_,q_2_, q_3_,v_1_,v_2_, v_3_);
   [generalized_force_term_] = GetGeneralizedForceTerm(w_H_1h_, dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_, q_1_, q_2_, q_3_, pressure_total(2,i), pressure_total(7,i), pressure_total(9,i));
   [kinetic_energy_term_] = GetKineticEnergyTerm(dw_H_1h__dq1_,dw_H_1h__dq2_,dw_H_1h__dq3_, w_H_2dot_);
   [potential_energy_term_] = GetPotentialEnergyTerm(q_1_, q_2_, q_3_, dw_H_1h__dq1_,dw_H_1h__dq2_, dw_H_1h__dq3_);
   
    
end