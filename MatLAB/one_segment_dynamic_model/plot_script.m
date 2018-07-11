
for i=1:numOfIterations
   time(i) = dt*i - dt;
   l_bar(i) = (q_(1,i) + q_(2,i) + q_(3,i))/3; 
   theta(i) = 2/3*(sqrt(q_(1,i)^2 + q_(2,i)^2 + q_(3,i)^2 - q_(1,i)*q_(2,i) - q_(1,i)*q_(3,i) - q_(2,i)*q_(3,i)))/0.026;
   
end
%%

figure; 
plot(time, q_(1,1:numOfIterations), time, q_(2,1:numOfIterations), time, q_(3,1:numOfIterations))
title('calculated muscle lengths p_{2} = 0.125,0.25,0.375,0.5,0')
xlabel('Time [s]')
ylabel('length [m]')

axis([0 50 0.25 0.35])
legend('q_{1}','q_{2}','q_{3}')

hold on 
plot(time, l_bar)
%%

figure;
plot(time, v_(1,1:numOfIterations), time, v_(2,1:numOfIterations) , time, v_(3,1:numOfIterations))
legend('v_1','v_2','v_3')
axis([0 25 0. 0.08])
xlabel('Time [s]')
ylabel('velocity [m/s]')
%%
figure;
plot(time,theta(1:90000))
xlabel('Time [s]')
ylabel('theta [rad]')
axis([0 25 0. 2])
%%


%%
figure;
plot(time, q_(2,1:numOfIterations))

%%
for i=1:50
    p(1,i) = 0;
    p(2,i) = 0;
    p(3,i) = i*0.5/50;
    t(i) = i/10;
end
for i=51:100
    p(1,i) = 0;
    p(2,i) = 0;
    p(3,i) = 0.5*(1-1/50*(i-50));
      t(i) = i/10;
end
for i=101:150
    p(1,i) = (i-100)*0.49/50;
    p(2,i) = 0;
    p(3,i) = (i-100)*0.5/50;
      t(i) = i/10;
end
for i=151:200
    p(1,i) = 0.49*(1-1/50*(i-150));
    p(2,i) = 0;
    p(3,i) = 0.5*(1-1/50*(i-150));
      t(i) = i/10;
end
%%
plot(t, p)
title('pressure input')
xlabel('time [s]')
ylabel('pressure [bar]')
legend('p_1', 'p_2', 'p_3')
axis([0 20 0 0.5])
%%
for i=1:50
    p(1,i) = 0;
    p(2,i) = 0.125;
    p(3,i) = 0;
    t(i) = i/10;
end
for i=51:100
    p(1,i) = 0;
    p(2,i) = 0.25;
    p(3,i) = 0;
      t(i) = i/10;
end
for i=101:150
    p(1,i) = 0;
    p(2,i) = 0.375;
    p(3,i) = 0;
      t(i) = i/10;
end
for i=151:200
    p(1,i) = 0;
    p(2,i) = 0.5;
    p(3,i) = 0;
      t(i) = i/10;
end
for i=201:250
    p(1,i) = 0;
    p(2,i) = 0;
    p(3,i) = 0;
      t(i) = i/10;
end

plot(t, p)
title('pressure input')
xlabel('time [s]')
ylabel('pressure [bar]')
legend('p_1', 'p_2', 'p_3')
axis([0 25 0 0.6])

