%% Cut and interpolate data
for i = 1:12000
    length_total_cut(:,i) = length_total(:,i);
    length_time_cut(i) = length_time(i);
end

length_time_cut_converted = seconds(length_time_cut);
pressure_time_converted = seconds(pressure_time);

length_time_table = timetable(length_time_cut_converted.', length_total_cut(1,:).', length_total_cut(2,:).', length_total_cut(3,:).');
length_time_table.Properties.VariableNames = {'s_1','s_2','s_3'};

pressure_time_table = timetable(pressure_time_converted, pressure_total(1,:).', pressure_total(2,:).', pressure_total(3,:).');
pressure_time_table.Properties.VariableNames = {'p1','p3','p2'};

interpolated_time_table = synchronize(pressure_time_table, length_time_table, 'regular', 'linear', 'SamplingRate', 100);
interpolated_time_table_fast = synchronize(pressure_time_table, length_time_table, 'regular', 'linear', 'SamplingRate', 200);

%% Filter static points
length_normalized(1,:) = interpolated_time_table_fast.s_1/max(interpolated_time_table_fast.s_1);
length_normalized(2,:) = interpolated_time_table_fast.s_2/max(interpolated_time_table_fast.s_2);
length_normalized(3,:) = interpolated_time_table_fast.s_3/max(interpolated_time_table_fast.s_3);

pressure_normalized(1,:) = interpolated_time_table_fast.p1/max(interpolated_time_table_fast.p1);

length_static_average = zeros(3,10);

for i = 0:9
    for j = 1:size(interpolated_time_table_fast,1)
        first_condition = (interpolated_time_table_fast.pressure_time_converted(j) > seconds((27.9 + i*10)));
        second_condition = (interpolated_time_table_fast.pressure_time_converted(j) < seconds((22.795 + (i+1)*10)));
        
        if(first_condition && second_condition)
            length_static_average(1, i+1) = length_static_average(1, i+1) + interpolated_time_table_fast.s_1(j);
            length_static_average(2, i+1) = length_static_average(2, i+1) + interpolated_time_table_fast.s_2(j);
            length_static_average(3, i+1) = length_static_average(3, i+1) + interpolated_time_table_fast.s_3(j);
            
            each_number_of_static_points(i+1) = j - last_else_j;
        else
            last_else_j = j; 
        end
    end
    length_static_average(:,i+1) = length_static_average(:,i+1)/each_number_of_static_points(i+1);  
end

for i=0:9
    pressure_static(i + 1) = (0.1 + i*0.05)*10^5;
end
%% main routine
%syms C k_spline;
%s_container = [441 401 385 356 337; 462 571 657 751 825; 435 413 374 347 320];
%p_container = [0.125; 0.25; 0.375; 0.5];
initial_length = [interpolated_time_table_fast.s_1(1); interpolated_time_table_fast.s_2(1); interpolated_time_table_fast.s_3(1)];

for i = 1:size(pressure_static,2)
    for j = 1:size(length_static_average,1)
        s_effective(i,j) = 0.3014 + (0.0384/1024)*(length_static_average(j,i) - initial_length(j));
    end
    
    [q_(i,:), theta_draw_wire_(i), phi_draw_wire(i), kappa_draw_wire(i)] = GetBalgLength(s_effective(i,:));
    
    [w_H_h_, theta_(i), phi_(i)] = GetTransformationMatrix(q_(i,:), theta_draw_wire_(i), phi_draw_wire(i), kappa_draw_wire(i));
    
    p_ = 10^5*[pressure_static(i); 0; 0];
    [parameters_(:,i), parameters_2_(:,i)] = SetEquations(theta_draw_wire_(i), phi_draw_wire(i)+pi/3, w_H_h_, q_, p_);
end

for i = 1:size(pressure_static,2)
    %C_set(i,1) = parameters_(i).C;
    %M_bending_set(i,1) = parameters_2_(i).M_bending_s;

    %k_spline_set(i,1) = parameters_(i).k_spline;
end
%%
theta_draw_wire_virtual = zeros(size(theta_draw_wire_, 2)+1, 1);
theta_draw_wire_virtual(1) = 0;

M_bending_set_virtual = zeros(size(theta_draw_wire_, 2)+1, 1);
M_bending_set_virtual(1) = 0;
for i = 2:(size(theta_draw_wire_, 2) + 1)
  theta_draw_wire_virtual(i)  = theta_draw_wire_(i-1);
  M_bending_set_virtual(i) = M_bending_set(i-1);
end
%%
plot(x_theta_phi_static(2,:), M_bending_set,'*')
% figure;
% plot(x_theta_phi_static(3,:),pressure_static)
% figure;
% plot(pressure_static, q,'*')
% axis( [1*10^4 6*10^4 0.27 0.35])

%% Define functions
function [parameters, parameters_2] = SetEquations(theta, phi, w_H_h, q, p)
    syms C M_bending_s;

    A = 0.0183*0.0183*pi; 
    k = 1108.63;
    r_b = 0.026; 
    m = 0.256; 
    g = 9.81; 
    q_0 = 0.3014;
    l_bar = (q(1) + q(2) + q(3))/3;
    
    h_r_B1 = [r_b; 0; 0; 1];
    h_r_B2 = [-r_b/2; (sqrt(3)/2)*r_b; 0; 1];
    h_r_B3 = [-r_b/2; -(sqrt(3)/2)*r_b; 0; 1];
    h_r_M  = [0; 0; 0; 1];
    
    w_r_B1 = w_H_h*h_r_B1;
    w_r_B2 = w_H_h*h_r_B2;
    w_r_B3 = w_H_h*h_r_B3;
    w_r_M = w_H_h*h_r_M;
    w_g = [0; 0; -g; 0];
    
    h_F_active_B1 = p(1)*A*[0; 0; 1; 0];
    h_F_active_B2 = p(2)*A*[0; 0; 1; 0];
    h_F_active_B3 = p(3)*A*[0; 0; 1; 0];

    h_F_passive_B1 = -k*(q(1) - q_0)*[0; 0; 1; 0];  
    h_F_passive_B2 = -k*(q(2) - q_0)*[0; 0; 1; 0];
    h_F_passive_B3 = -k*(q(3) - q_0)*[0; 0; 1; 0];

    h_F_B1 = h_F_active_B1 + h_F_passive_B1; 
    h_F_B2 = h_F_active_B2 + h_F_passive_B2;
    h_F_B3 = h_F_active_B3 + h_F_passive_B3;
    
    %h_F_spline = -k_spline*(l_bar - q_0)*[0; 0; 1; 0];
    
    w_F_B1 = w_H_h*h_F_B1;
    w_F_B2 = w_H_h*h_F_B2;
    w_F_B3 = w_H_h*h_F_B3;
    %w_F_spline = w_H_h*h_F_spline;
    
    M_forces = cross(w_F_B1(1:3), w_r_B1(1:3)) + cross(w_F_B2(1:3), w_r_B2(1:3)) + cross(w_F_B3(1:3), w_r_B3(1:3)); %+ cross(w_F_spline(1:3), w_r_M(1:3));
    M_gravity = cross(m*w_g(1:3), w_r_M(1:3));
    M_bending = -C*theta*[-sin(phi); cos(phi); 0];
    M_bending_2 = M_bending_s*[-sin(phi); cos(phi); 0];
    
    eqns = M_forces + M_gravity + M_bending == [0; 0; 0];
    eqns_2 = M_forces + M_gravity + M_bending_2 == [0; 0; 0];
    [A,b] = equationsToMatrix(eqns(1:2), C)
    [A_2,b_2] = equationsToMatrix(eqns_2(1:2), M_bending_s)
    
    parameters = mldivide(A,b)
    parameters_2 = A_2\b_2
end

function [q, theta, phi, kappa] = GetBalgLength(s)
    phi_0 = pi/3;
    r_b = 0.026; 
    r_s = 0.010;
    
    theta = (2/3)*((sqrt(s(1)^2 + s(2)^2 + s(3)^2 - s(1)*s(2) - s(1)*s(3) - s(2)*s(3)))/r_s);
    phi = atan2((sqrt(3)*(s(3) - s(2))),(s(2) + s(3) - 2*s(1)));
    l_bar = (s(1) + s(2) + s(3))/3;
    kappa = theta/l_bar;
   
    if theta == 0
        q(1) = 0.3014;
        q(2) = 0.3014;
        q(3) = 0.3014;
    else
        q(1) = theta*(1/kappa - r_b*cos(phi + phi_0));
        q(2) = theta*(1/kappa - r_b*cos(phi + 2*pi/3 + phi_0));
        q(3) = theta*(1/kappa - r_b*cos(phi + 4*pi/3 + phi_0));
    end
end

function [w_H_h, phi, theta] = GetTransformationMatrix(q, theta, phi, kappa)
    

    w_H_h_1 = [cos(phi)^2*(cos(theta) - 1) + 1, sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)*sin(theta), (cos(phi)*(1 - cos(theta)))/kappa];
    w_H_h_2 = [sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)^2*(1 - cos(theta)) + cos(theta), sin(phi)*sin(theta), (sin(phi)*(1 - cos(theta)))/kappa];
    w_H_h_3 = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta), sin(theta)/kappa];
    w_H_h_4 = [0, 0, 0, 1];
    w_H_h = [w_H_h_1; w_H_h_2; w_H_h_3; w_H_h_4];
end
