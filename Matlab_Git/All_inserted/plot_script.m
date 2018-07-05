figure; 
plot(time, q_(1,1:numOfIterations), time, q_(2,1:numOfIterations), time, q_(3,1:numOfIterations))
title('calculated muscle lengths p2 = 0.25')
xlabel('Time [s]')
ylabel('length [m]')
legend('q1','q2','q3')
axis([0 0.5 0 0.6])

