%This is a subroutine within the subroutine 'delta13C_soilCO2' within PBUQ

%Compiled data from Quade et al (in review) and Passey et al (2010)

D47T = [28.40625; 22.052; 26.375; 35.5; 32; 30; 33; 30.7; 24.175; 21.2; 20.175; 17.895; 11.23; 23.67; 19.76; 26; 32.67142857; 18.8; 23.125; 26.06666667];
Mean_annual_air_temperature = [18.2; 11; 6.7; 22.9; 25.6; 25.5; 24.4; 23.5; 7.6; 7.6; 3; 5; -1.6; 3.6; 4.6; 27.6; 30.3; 16.5; 6; 15.3];

yhat = 0.506 * Mean_annual_air_temperature + 17.974;

SSE = sum((yhat-D47T).^2);

lengthD47T = length(D47T);

MSE = SSE/(lengthD47T-2);

mean_mean_annual_air_temperature = mean(Mean_annual_air_temperature);

sum_xi_minus_xbar_squared = sum((Mean_annual_air_temperature-mean_mean_annual_air_temperature).^2);

all_MAAT = linspace (-5,35,201);

StandardError = sqrt(MSE*(1/lengthD47T + (all_MAAT-mean_mean_annual_air_temperature).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 202

StandardDeviation = StandardError / (sqrt(1/lengthD47T));

StandardError_newobs = sqrt(MSE*(1 + 1/lengthD47T + (all_MAAT-mean_mean_annual_air_temperature).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 204
%the variance of new indivual observation will be larger than the variance
%of means of observations (the variance of the means is related to the 
%standard error above). The variance in a new observation is what we are
%interested in

yhat_all = 0.506 * all_MAAT + 17.974;

%Below plots data and regression line and confidence intervals
% figure(7)  
% plot (Mean_annual_air_temperature, D47T, 'ko')
% hold on
% plot (all_MAAT, yhat_all)
% plot (all_MAAT, yhat_all + StandardError, 'b--')
% plot (all_MAAT, yhat_all - StandardError, 'b--')
% plot (all_MAAT, yhat_all + StandardDeviation, 'b--')
% plot (all_MAAT, yhat_all - StandardDeviation, 'b--')
% plot (all_MAAT, yhat_all + StandardError_newobs, 'r--')
% plot (all_MAAT, yhat_all - StandardError_newobs, 'r--')
% xlabel('MAAT')
% ylabel('D47 T')