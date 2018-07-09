%% Define variables
syms q_1_dot q_2_dot q_3_dot v_1_dot v_2_dot v_3_dot;

%% Parameter list. !!After Changing parameters, run this section.!! 
D_damp_spline_ = 100;
m_ = 0.265;
g_ = -9.81;
C_ohne_l_bar_ = -0.26;
q_0_ = 0.3014;
A_ = 0.0183*0.0183*pi;
k_ =    616.95;
D_damp_ = 200;
r_b_ = 0.026;
k_spline_ = 2.9822*10^4;

%% main section
numOfIterations = 90000 ;
dt = 0.0005;

q_ = zeros(3, numOfIterations);
v_ = zeros(3, numOfIterations);
p_ = zeros(3, numOfIterations);
time = zeros(1, numOfIterations);

q_(:,1) = [0.3014; 0.3014; 0.3014];
v_(:,1) = [0; 0; 0];
for i=1:1000 
p_(1,i) = 0.*10^5;
p_(2,i) = 0.25*10^5;  
p_(3,i) = 0.*10^5; 
 end
for i=10001 : 20000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.25*10^5;  
p_(3,i) = 0.*10^5; 
end
for i=20001 : 30000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.375*10^5;  
p_(3,i) = 0.*10^5; 
end
for i=30001 : 40000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.5*10^5;  
p_(3,i) = 0.*10^5; 
end

for i=40001 : 50000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.*10^5;  
p_(3,i) = 0.*10^5; 
end

for i=50001 : 60000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.125*10^5;  
p_(3,i) = 0.125*10^5; 
end
for i=60001 : 70000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.25*10^5;  
p_(3,i) = 0.25*10^5; 
end
for i=70001 : 80000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.375*10^5;  
p_(3,i) = 0.375*10^5; 
end

for i=80001 : 90000
p_(1,i) = 0.*10^5;
p_(2,i) = 0.5*10^5;  
p_(3,i) = 0.5*10^5; 
end



% for i=20100:30000
%     p_(1,:) = 0.*10^5;
% p_(2,:) = 0.*10^5;  
% p_(3,:) = 0.25*10^5; 
% end
% for i=30100:40000
%     p_(1,:) = 0.25*10^5;
% p_(2,:) = 0.25*10^5;  
% p_(3,:) = 0.25*10^5; 
% end
% 
% for i=30100:40000
%     p_(1,:) = 0.25*10^5;
% p_(2,:) = 0.25*10^5;  
% p_(3,:) = 0.*10^5; 
% end
%p_(:,:) = [0; 0.25*10^5; 0.25*10^5];
% for i=40100:50000
%     p_(1,:) = 0.5*10^5/10000*(i-401000);
% p_(2,:) = 0.5*10^5/10000*(i-401000);  
% p_(3,:) = 0.*10^5; 
% end
% for i=50100:60000
%     p_(1,:) = 0.*10^5/10000*(i-401000);
% p_(2,:) = 0.*10^5/10000*(i-401000);  
% p_(3,:) = 0.*10^5; 
% end
% 
% for i=60100:70000
% p_(1,:) = 0.*10^5/10000*(i-401000);
% p_(2,:) = 0.*10^5/10000*(i-401000);  
% p_(3,:) = 0.6*10^5; 
% end
% 
% for i=70100:80000
% p_(1,:) = 0.*10^5/10000*(i-401000);
% p_(2,:) = 0.2*10^5/10000*(i-401000);  
% p_(3,:) = 0.6*10^5; 
% end
% for i=80100:100000
% p_(1,:) = 0.*10^5/10000*(i-401000);
% p_(2,:) = 0.6 *10^5/10000*(i-401000);  
% p_(3,:) = 0.*10^5; 
% end
k_factors_ = zeros(4,6);

for i = 1:numOfIterations
    k_factors_ = zeros(4,6);
    
    [k_factors_(1,1:3), k_factors_(1,4:6)] = (SolveSystemRungeKutta(D_damp_spline_, m_, g_, C_ohne_l_bar_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i), q_(2,i), q_(3,i), v_(1,i), v_(2,i), v_(3,i), q_1_dot, q_2_dot, q_3_dot, v_1_dot, v_2_dot, v_3_dot,p_(1,i), p_(2,i), p_(3,i)));
    k_factors_(1,:) = dt*k_factors_(1,:);
    
    [k_factors_(2,1:3), k_factors_(2,4:6)] = (SolveSystemRungeKutta(D_damp_spline_, m_, g_, C_ohne_l_bar_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i) + k_factors_(1,1)/2, q_(2,i) + k_factors_(1,2)/2, q_(3,i) + k_factors_(1,3)/2, v_(1,i) + k_factors_(1,4)/2, v_(2,i) + k_factors_(1,5)/2, v_(3,i) + k_factors_(1,6)/2, q_1_dot, q_2_dot, q_3_dot, v_1_dot, v_2_dot, v_3_dot, p_(1,i), p_(2,i), p_(3,i)));
    k_factors_(2,:) = dt*k_factors_(2,:);
    
    [k_factors_(3,1:3), k_factors_(3,4:6)] = (SolveSystemRungeKutta(D_damp_spline_, m_, g_, C_ohne_l_bar_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i) + k_factors_(2,1)/2, q_(2,i) + k_factors_(2,2)/2, q_(3,i) + k_factors_(2,3)/2, v_(1,i) + k_factors_(2,4)/2, v_(2,i) + k_factors_(2,5)/2, v_(3,i) + k_factors_(2,6)/2, q_1_dot, q_2_dot, q_3_dot, v_1_dot, v_2_dot, v_3_dot, p_(1,i), p_(2,i), p_(3,i)));
    k_factors_(3,:) = dt*k_factors_(3,:);
    
    [k_factors_(4,1:3), k_factors_(4,4:6)] = (SolveSystemRungeKutta(D_damp_spline_, m_, g_, C_ohne_l_bar_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i) + k_factors_(3,1), q_(2,i) + k_factors_(3,2), q_(3,i) + k_factors_(3,3), v_(1,i) + k_factors_(3,4), v_(2,i) + k_factors_(3,5), v_(3,i) + k_factors_(3,6), q_1_dot, q_2_dot, q_3_dot, v_1_dot, v_2_dot, v_3_dot, p_(1,i), p_(2,i), p_(3,i)));
    k_factors_(4,:) = dt*k_factors_(4,:);
    
    q_(:,i+1) = q_(:,i) + 1/6*(k_factors_(1,1:3).' + 2*k_factors_(2,1:3).' + 2*k_factors_(3,1:3).' + k_factors_(4,1:3).');
    v_(:,i+1) = v_(:,i) + 1/6*(k_factors_(1,4:6).' + 2*k_factors_(2,4:6).' + 2*k_factors_(3,4:6).' + k_factors_(4,4:6).');
end
