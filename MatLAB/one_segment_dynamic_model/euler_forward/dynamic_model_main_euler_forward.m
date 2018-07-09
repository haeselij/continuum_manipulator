%% Define variables
syms q_1_n q_2_n q_3_n v_1_n v_2_n v_3_n
assume([q_1_n, q_2_n, q_3_n, v_1_n, v_2_n, v_3_n], 'real');
%% Parameter list. !!After Changing parameters, run this section.!! 
D_damp_spline_ = 100;
m_ = 0.265;
g_ = -9.81;
C_ohne_l_bar_ = -0.15;
q_0_ = 0.3014;
A_ = 0.0183*0.0183*pi;
k_ = 616.95;
D_damp_ = 200;
r_b_ = 0.026;
k_spline_ = 2.9822*10^4;

%% main section
numOfIterations = 20000;
dt = 0.0005;

q_ = zeros(3, numOfIterations);
v_ = zeros(3, numOfIterations);
p_ = zeros(3, numOfIterations);
time = zeros(1, numOfIterations);

q_(:,1) = [0.3014; 0.3014; 0.3014];
v_(:,1) = [0; 0; 0];
for i=1:5000
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5/5000*i;
    p_(3,i) = 0;
end

for i=5001:10000
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5*(1 - 1/5000*(i-5000));
    p_(3,i) = 0; 
end
for i=10001:15000
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5/5000*(i-10000);
    p_(3,i) = 0.5*10^5/5000*(i-10000);
end

for i=15001:20000
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5*(1 - 1/5000*(i-15000));
    p_(3,i) = 0.5*10^5*(1 - 1/5000*(i-15000));
end

for i = 1:numOfIterations
    tic;
    q_1_dot = 1/dt*(q_1_n - q_(1,i));
    q_2_dot = 1/dt*(q_2_n - q_(2,i));
    q_3_dot = 1/dt*(q_3_n - q_(3,i));
    
    v_1_dot = 1/dt*(v_1_n - v_(1,i));
    v_2_dot = 1/dt*(v_2_n - v_(2,i));
    v_3_dot = 1/dt*(v_3_n - v_(3,i));
    
    [q_(:,i+1), v_(:,i+1)] = SolveSystemEulerForward(dt ,D_damp_spline_, m_, g_, C_ohne_l_bar_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i), q_(2,i), q_(3,i), v_(1,i), v_(2,i), v_(3,i), q_1_dot, q_2_dot, q_3_dot, v_1_dot, v_2_dot, v_3_dot, p_(1,i), p_(2,i), p_(3,i));
    time(1,i) = dt*i - dt;
    toc
end


