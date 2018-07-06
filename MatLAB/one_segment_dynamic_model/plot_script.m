for i=1:numOfIterations
   time(i) = dt*i - dt;
   l_bar(i) = (q_(1,i) + q_(2,i) + q_(3,i))/3; 
   theta(i) = 2/3*(sqrt(q_(1,i)^2 + q_(2,i)^2 + q_(3,i)^2 - q_(1,i)*q_(2,i) - q_(1,i)*q_(3,i) - q_(2,i)*q_(3,i)))/0.026;
   
end
%%
figure; 
plot(time, q_(1,1:numOfIterations), time, q_(2,1:numOfIterations), time, q_(3,1:numOfIterations))
title('calculated muscle lengths p_{2} = 0.5')
xlabel('Time [s]')
ylabel('length [m]')

axis([0 4.5 0.25 0.35])
legend('q_{1}','q_{2}','q_{3}')
%%
hold on 
plot(time, l_bar)
%%

figure;
plot(time, v_(1,1:numOfIterations), time, v_(2,1:numOfIterations) , time, v_(3,1:numOfIterations))
legend('v1','v2','v3')
%%
figure;
plot(time,theta(1:10000))
xlabel('Time [s]')
ylabel('theta [rad]')
%%


%%
figure;
plot(time, q_(2,1:numOfIterations))
