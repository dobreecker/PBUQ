%data from Cotton&Sheldon(in press)

%depth_to_Bk_Retallack_2009 = [110, 78, 66, 40, 40, 40, 30, 36, 40, 28, 40, 20, 13, 13, 35];
%Sz_Retallack_2009 = [7921, 5616, 4342, 4474, 4274, 4614, 2374, 2853, 2274, 2417, 3938, 2339, 1404, 664, 1153];

%revised to reflect minimum warm season S(z)
depth_to_Bk_Retallack_2009 = [110, 78, 66, 50, 50, 50, 30, 36, 40, 28, 40, 13, 13, 35, 20, 2, 2];
Sz_Retallack_2009 = [4185, 3085, 1876, 3172, 3375, 4172, 1471, 1800, 1876, 957, 1685, 1401, 751, 1063, 1591, 611, 325];

%yhat = 67.0*depth_to_Bk_Retallack_2009+567;

%revised to reflect minimum warm season S(z)
yhat = 35.3*depth_to_Bk_Retallack_2009+588;

SSE = sum((yhat-Sz_Retallack_2009).^2);

length_Sz_Retallack_2009 = length(Sz_Retallack_2009);

MSE = SSE/(length_Sz_Retallack_2009-2);

mean_Retallack_2009_depth_to_Bk = mean(depth_to_Bk_Retallack_2009);

sum_xi_minus_xbar_squared = sum((depth_to_Bk_Retallack_2009-mean_Retallack_2009_depth_to_Bk).^2);

all_depth_to_Bk = linspace(0,120,121);

StandardError = sqrt(MSE*(1/length_Sz_Retallack_2009 + (all_depth_to_Bk-mean_Retallack_2009_depth_to_Bk).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 202

StandardDeviation = StandardError / (sqrt(1/length_Sz_Retallack_2009));

StandardError_newobs = sqrt(MSE*(1 + 1/length_Sz_Retallack_2009 + (all_depth_to_Bk-mean_Retallack_2009_depth_to_Bk).^2/sum_xi_minus_xbar_squared));
%the variance of new indivual observation will be larger than the variance
%of means of observations (the variance of the means is related to the 
%standard error above). The variance in a new observation is what we are
%interested in


%yhat_all = 67.0*all_depth_to_Bk+567;

%revised to reflect minimum warm season S(z)
yhat_all = 35.3*all_depth_to_Bk+588;


%Below plots data and regression line and confidence intervals
% figure(6)
% plot (depth_to_Bk_Retallack_2009, Sz_Retallack_2009, 'ko')
% hold on
% plot (all_depth_to_Bk, yhat_all)
% % plot (all_depth_to_Bk, yhat_all + StandardError, 'b--')
% % plot (all_depth_to_Bk, yhat_all - StandardError, 'b--')
% % plot (all_depth_to_Bk, yhat_all + StandardDeviation, 'b--')
% % plot (all_depth_to_Bk, yhat_all - StandardDeviation, 'b--')
% plot (all_depth_to_Bk, yhat_all + StandardError_newobs, 'r--')
% plot (all_depth_to_Bk, yhat_all - StandardError_newobs, 'r--')
% xlabel('Depth to Bk (cm)')
% ylabel('S(z) (ppmV)')



%below plots the relative error versus absolute value of S(z)
%figure(7)
%plot (yhat_all, StandardError_newobs./yhat_all)
