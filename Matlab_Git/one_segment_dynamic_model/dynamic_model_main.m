p_container = zeros(3, pressures_size);  
q_container= zeros(3, pressures_size);  


w_H_1h = zeros(4,4); % coordinatetransform matrix head to worldframe

d_w_H_1h__d_q_1 = zeros(4,4); %first derrivative of H_1h after seg.q_1
d_w_H_1h__d_q_2 = zeros(4,4);
d_w_H_1h__d_q_3 = zeros(4,4);

w_H_1h_2dot = zeros(4,1);