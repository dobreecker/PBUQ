%This is a subroutine within the subroutine 'Monte_Carlo_error_propagation' within PBUQ

load Kohn_2010_data.txt

log_MAP_plus_300 = Kohn_2010_data(:,1);
BIG_DELTA = Kohn_2010_data(:,2);

yhat = 5.9634 * log_MAP_plus_300 + 1.9539;

SSE = sum((yhat-BIG_DELTA).^2);

length_BIG_DELTA = length(BIG_DELTA);

MSE = SSE/(length_BIG_DELTA-2);

mean_log_MAP_plus_300 = mean(log_MAP_plus_300);

sum_xi_minus_xbar_squared = sum((log_MAP_plus_300-mean_log_MAP_plus_300).^2);

all_log_MAP_plus_300 = linspace(2,4,100);

StandardError = sqrt(MSE*(1/length_BIG_DELTA + (all_log_MAP_plus_300-mean_log_MAP_plus_300).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 202


StandardDeviation = StandardError / (sqrt(1/length_BIG_DELTA));

StandardError_newobs = sqrt(MSE*(1 + 1/length_BIG_DELTA + (all_log_MAP_plus_300-mean_log_MAP_plus_300).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 204
%the variance of new indivual observation will be larger than the variance
%of means of observations (the variance of the means is related to the 
%standard error above). The variance in a new observation is what we are
%interested in

yhat_all = 5.9634 * all_log_MAP_plus_300 + 1.9539;

%Below plots data and regression line and confidence intervals
% figure(8)
% plot (log_MAP_plus_300, BIG_DELTA, 'ko')
% hold on
% plot (all_log_MAP_plus_300, yhat_all)
% plot (all_log_MAP_plus_300, yhat_all + StandardError, 'b--')
% plot (all_log_MAP_plus_300, yhat_all - StandardError, 'b--')
% plot (all_log_MAP_plus_300, yhat_all + StandardDeviation, 'b--')
% plot (all_log_MAP_plus_300, yhat_all - StandardDeviation, 'b--')
% plot (all_log_MAP_plus_300, yhat_all + StandardError_newobs, 'r--')
% plot (all_log_MAP_plus_300, yhat_all - StandardError_newobs, 'r--')
% xlabel('log (MAP+300)')
% ylabel('Big Delta (atm-leaf)')