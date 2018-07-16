
%% Parameter list. !!After Changing parameters, run this section.!! 
D_damp_spline_ = 100;
m_ = 0.256;
g_ = -9.81;
C_ = -0.15;
q_0_ = 0.3014;
A_ = 0.0183*0.0183*pi;
k_ = 616.95;

D_damp_ = 200;
r_b_ = 0.026;
k_spline_ = 2.9822*10^4;

%% main section
dt = 0.0005;
numOfIterations = 20/dt;
numDataPer5Sek = 5/dt;

q_ = zeros(3, numOfIterations);
v_ = zeros(3, numOfIterations);
p_ = zeros(3, numOfIterations);
time = zeros(1, numOfIterations);

q_(:,1) = [0.3014; 0.3014; 0.3014];
v_(:,1) = [0; 0; 0];

for i=1:numDataPer5Sek
    
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5/numDataPer5Sek*i;
    p_(3,i) = 0;
    
end

for i= numDataPer5Sek+1:2*numDataPer5Sek
    
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5*(1 - 1/numDataPer5Sek*(i-numDataPer5Sek));
    p_(3,i) = 0; 
    
end

for i=2*numDataPer5Sek+1:3*numDataPer5Sek
    
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5/numDataPer5Sek*(i-2*numDataPer5Sek);
    p_(3,i) = 0.5*10^5/numDataPer5Sek*(i-2*numDataPer5Sek);
    
end

for i=3*numDataPer5Sek+1:4*numDataPer5Sek
    
    p_(1,i) = 0;
    p_(2,i) = 0.5*10^5*(1 - 1/numDataPer5Sek*(i-3*numDataPer5Sek));
    p_(3,i) = 0.5*10^5*(1 - 1/numDataPer5Sek*(i-3*numDataPer5Sek));
    
end

tic;
for i = 1:numOfIterations-1
    
    [q_(:,i+1), v_(:,i+1)] = SolveSystemEulerForward(dt, D_damp_spline_, m_, g_, C_, q_0_, A_, k_, D_damp_, r_b_, k_spline_, q_(1,i), q_(2,i), q_(3,i), v_(1,i), v_(2,i), v_(3,i), p_(1,i), p_(2,i), p_(3,i));
    time(1,i) = dt*i - dt;
    
end
toc

