
%% Cut and read static points 
j = 1;
k = 1;
time = 0;

for i= 1:11630
    
   if length_time(i) > 13.5 + j
       
    length(1,k) = (length_total(1,i) - length_total(1,1))*(0.0384/1024) + 0.3015;
    length(2,k) = (length_total(2,i) - length_total(1,2))*(0.0384/1024) + 0.3015;
    length(3,k) = (length_total(3,i)- length_total(1,1))*(0.0384/1024) + 0.3015;   
    pressure_static(k) = (0.05 + (k-1)*0.05)*10^5;
    j = j + 10;
    k = k + 1;
     
   end
end

for i= 0:10
    
    pressure_static(i + 1) = (0.1 + i*0.05)*10^5;
    
end

x_theta_phi_static = zeros(3,11);
C_set = zeros(1,11);
M_bending_set = zeros(1,11);

for i = 1:11
    [x_theta_phi_static(1,i), x_theta_phi_static(2,i), x_theta_phi_static(3,i)] = GetXThetaPhi(length(1,i), length(2,i),length(3,i));
    [q_1, q_2, q_3] = GetBalgLength(x_theta_phi_static(2,i), x_theta_phi_static(3,i));
    q(i) = q_1;
    [C_set(i)] = GetBendingStiffnessC(q_1, q_2, q_3, 0.3014, 0.3014, 0.3014, pressure_static(i), 0, 0, x_theta_phi_static(1,i), x_theta_phi_static(2,i));
end

C_mean = mean(C_set);

plot(pressure_static*10 ^-5,x_theta_phi_static(2,:),'*')
xlabel('pressure [bar]')
ylabel('\theta [rad]')
set(gca,'fontsize', 12);

%% Define functions

function [x,theta, phi] = GetXThetaPhi(q_1, q_2, q_3)

    theta = (2/3)*((sqrt(q_1^2 + q_2^2 + q_3^2 - q_1*q_2 - q_1*q_3 - q_2*q_3))/0.01);
    phi = atan2((sqrt(3)*(q_3 - q_2)),(q_2 + q_3 - 2*q_1));
    l_bar = (q_1 + q_2 + q_3)/3.0;
    kappa = theta/l_bar;
    
    x = cos(phi)*(1-cos(theta))/kappa;
end

function [C] = GetBendingStiffnessC(q1, q2, q3, q10, q20, q30, p1, p2, p3, x, theta)

    syms C_unknown;

    ratio_A_k = 0.095/(10^5); 
    A = 0.0163*0.0163*pi; 
    k = 1108.8;
    r_b = 0.026; 
    m = 0.256; 
    g = 9.81;
    l_bar = 0.3014;
    kappa = theta/l_bar;

    equation1 = m*g*cos(theta)*1/kappa*(1-cos(theta/2)) + r_b*(-k*(q1 - q10) + p1*A) - 1/2*(-k*(q2 - q20) + p2*A)*r_b - 1/2*(-k*(q3 - q30) + p3*A)*r_b - C_unknown*10*theta == 0;

    C = solve(equation1, C_unknown);
  
end



function [q_11, q_12, q_13] = GetBalgLength(theta, phi)

    phi_0 = -2*pi/3;
    l_bar = 0.3014;
    r_b = 0.026;
    kappa = theta/l_bar;
   
    if theta == 0
        q_11 = 0.3014;
        q_12 = 0.3014;
        q_13 = 0.3014;
    else
        q_11 = theta*(1/kappa - r_b*cos(phi + phi_0));
        q_12 = theta*(1/kappa - r_b*cos(phi + 2*pi/3 + phi_0));
        q_13 = theta*(1/kappa - r_b*cos(phi + 4*pi/3 + phi_0));
        
    end
end
