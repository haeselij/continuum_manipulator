%%
for i=1:numOfIterations
   time(i) = 2*dt*i - dt; 
   l_bar_container(i) = (q_(1,i) + q_(2,i) + q_(3,i))/3;
   
end

figure; 
plot(time, q_(3,1:numOfIterations), time, q_(1,1:numOfIterations), time, q_(2,1:numOfIterations))
title('calculated muscle lengths euler  ')
xlabel('Time [s]')
ylabel('Length [m]')
legend('q_1','q_2','q_3')
axis([0 20 0.25 0.35])
hold on
%plot(time, l_bar_container)
%%
figure;
e_r = q_euler - q_;

plot(time,e_r(1:20000))