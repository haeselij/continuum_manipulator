%% Filter out steady state values

j = 1;
t = 0;
t_2 = 0;
p = 1;

length_total_delta = length_total*(0.0384/1024) - 0.015;

    for i=1:length_size
        
        if(length_time(i) > t + 4.7 )
           length_stiffness(j)  = length_total_delta(1,i);
           j = j+ 1;
           t= t + 10;
           
        end

    end
    
    for i =1: pressure_size

        if (pressure_time(i) > t_2 + 4.7)
            pressure_stiffness(p) = pressure_total(2,i)*10^5;
            p = p + 1;
            t_2 = t_2 + 10;
            
        end
        
    end

figure;
plot(pressure_stiffness(2:12), length_stiffness(2:12),'*')
title('longitudioal stiffness measurement')
xlabel('pressure [bar]')
ylabel(' \delta draw wire measurement [mm]')

%% calculate k

p_l_ratio = polyfit(pressure_stiffness(2:12),length_stiffness(2:12),1);

k_M = (pi * 0.0183^2)/(p_l_ratio(1));


