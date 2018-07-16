%% Analize data

%%  interpolate data
exp_1_start = 1070;
exp_1_end = 4722;

exp_2_start = 6025; %-32.82
exp_2_end = 9596;

exp_3_start = 14353;
exp_3_end = 18040;
q_exp_1(:,:) = q_exp(:,exp_1_start:exp_1_end);
q_model_1(:,:) = q_; 

exp_1_time = length_time(exp_1_start:exp_1_end) - length_time(exp_1_start);
exp_1_time_converted = seconds(exp_1_time);

for i = 1 : size(q_model_1, 2)
   model_1_time(i) = 0.001*i - 0.001; 
end
model_1_time_converted = seconds(model_1_time.');

ex_1_timetable = timetable(exp_1_time_converted, q_exp_1(1,:).',q_exp_1(2,:).', q_exp_1(3,:).');
mod_1_timetable = timetable(model_1_time_converted, q_model_1(1,:).',q_model_1(2,:).', q_model_1(3,:).');

ex_1_timetable.Properties.VariableNames = {'q_1exp','q_2exp','q_3exp'};
mod_1_timetable.Properties.VariableNames = {'q_1mod','q_2mod','q_3mod'};

sync_time_table = synchronize(ex_1_timetable,mod_1_timetable, 'regular', 'linear', 'SamplingRate', 100);
muscle_error_vector(1,:) = sync_time_table.q_1exp - sync_time_table.q_1mod;
muscle_error_vector(2,:) = sync_time_table.q_2exp - sync_time_table.q_2mod;
muscle_error_vector(3,:) = sync_time_table.q_3exp - sync_time_table.q_3mod;


%%
size_sync_time_table = size(sync_time_table,1);
for i =1: size_sync_time_table
    
[exp_tip(1,i), exp_tip(2,i), exp_tip(3,i), exp_tip(4,i), exp_theta(i), exp_phi(i), exp_l_bar(i)] = get_tip(sync_time_table.q_1exp(i), sync_time_table.q_2exp(i), sync_time_table.q_3exp(i));
[mod_tip(1,i),mod_tip(2,i), mod_tip(3,i), mod_tip(4,i), mod_theta(i),mod_phi(i), mod_l_bar(i)] = get_tip(sync_time_table.q_1mod(i), sync_time_table.q_2mod(i),sync_time_table.q_3mod(i));

end 
coordinate_error_vector(:,:) = exp_tip - mod_tip;
for i= 1:size_sync_time_table
position_error_absv(i) = norm(coordinate_error_vector(1:3,i));
muscle_error_absv(i) = norm(muscle_error_vector(:,i));

end
%%
position_error_mean = sum(position_error_absv)/size_sync_time_table
muscle_error_mean = sum(muscle_error_absv)/size_sync_time_table
muscle_std = std(muscle_error_vector.')
position_std = std(position_error_absv')
%%
figure;
plot(sync_time_table.exp_1_time_converted, exp_tip(1,:),sync_time_table.exp_1_time_converted, mod_tip(1,:))
title('comparison x coordinate')
xlabel('time [s]')
ylabel('x [m]')
legend('experiment', 'model')
axis([0.0 20.0 -0.2  0.1])

%%
figure;
plot(sync_time_table.exp_1_time_converted, exp_tip(2,:),sync_time_table.exp_1_time_converted, mod_tip(2,:))
title('comparison y coordinate')
xlabel('time [s]')
ylabel('y [m]')
legend('experiment', 'model')
axis([0.0 20.0 -0.2  0.1])
%%
figure;
plot(sync_time_table.exp_1_time_converted, exp_tip(3,:),sync_time_table.exp_1_time_converted, mod_tip(3,:))
title('comparison z coordinate')
xlabel('time [s]')
ylabel('z [m]')
legend('experiment', 'model')

%%
% figure;
% plot(sync_time_table.exp_1_time_converted,muscle_error_vector)
% %%
figure;
plot(sync_time_table.exp_1_time_converted, sync_time_table.q_1exp, sync_time_table.exp_1_time_converted, sync_time_table.q_2exp, sync_time_table.exp_1_time_converted, sync_time_table.q_3exp)
%% theta comparison

figure;
plot(sync_time_table.exp_1_time_converted, exp_theta, sync_time_table.exp_1_time_converted, mod_theta)
title('comparison z coordinate')
xlabel('time [s]')
ylabel('z [m]')
legend('experiment', 'model')

%%
plot(sync_time_table.exp_1_time_converted, sync_time_table.q_1mod, sync_time_table.exp_1_time_converted, sync_time_table.q_2mod, sync_time_table.exp_1_time_converted, sync_time_table.q_3mod)

%%
function [x,y,z,w, theta, phi, l_bar] = get_tip(q_11, q_12, q_13) 
    if q_11 == q_12 && q_12 == q_13
        x = 0;
        y= 0;
        z = .3014;
        w= 1;
        
        theta = 0;
        phi = 0;
        l_bar = (q_11 + q_12 + q_13)/3;
    else 
        theta = 2/3*(sqrt(q_11^2 + q_12^2 + q_13^2 - q_11*q_12 - q_11*q_13 - q_12*q_13))/0.026;
        l_bar = (q_11 + q_12 + q_13)/3;
        phi = atan2((sqrt(3)*(q_13 - q_12)),(q_12 + q_13 - 2*q_11));
        kappa = theta/l_bar;

        w_H_1h_1 = [cos(phi)^2*(cos(theta) - 1) + 1, sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)*sin(theta), (cos(phi)*(1 - cos(theta)))/kappa];
        w_H_1h_2 = [sin(phi)*cos(phi)*(cos(theta) - 1), cos(phi)^2*(1 - cos(theta)) + cos(theta), sin(phi)*sin(theta), (sin(phi)*(1 - cos(theta)))/kappa];
        w_H_1h_3 = [-cos(phi)*sin(theta), -sin(phi)*sin(theta), cos(theta), sin(theta)/kappa];
        w_H_1h_4 = [0, 0, 0, 1];
        w_H_1h = [w_H_1h_1; w_H_1h_2; w_H_1h_3; w_H_1h_4];
        
        
        
        %w_H_t_total(:,:,i) = R*w_H_1h;
        r_tip_pos = w_H_1h*[0; 0; 0; 1];
 
        y = - r_tip_pos(1);
        x = r_tip_pos(2);
        z = r_tip_pos(3);
        w = r_tip_pos(4);
    end
end