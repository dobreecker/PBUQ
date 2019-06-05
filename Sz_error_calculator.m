%data from Cotton&Sheldon(in press)
Cotton_SheldonMAP = [153; 114; 350; 278; 375; 210; 220; 374; 330; 339; 301; 605; 158; 312; 116; 124; 225; 540; 484; 490; 345; 289; 495];
Cotton_SheldonSz = [554; 354; 1471; 876; 1988; 669; 692; 318; 3376; 1872; 1685; 3287; 1063; 892; 315; 425; 853; 2655; 3085; 2174; 3172; 957; 2059];

%meanCandSSz = mean(Cotton_SheldonSz);

yhat = 5.673*Cotton_SheldonMAP-269.87;

SSE = sum((yhat-Cotton_SheldonSz).^2);

length_Cotton_SheldonSz = length(Cotton_SheldonSz);

MSE = SSE/(length_Cotton_SheldonSz-2);

meanCandSMAP = mean(Cotton_SheldonMAP);

sum_xi_minus_xbar_squared = sum((Cotton_SheldonMAP-meanCandSMAP).^2);

all_MAP = linspace(0,700,701);

StandardError = sqrt(MSE*(1/length_Cotton_SheldonSz + (all_MAP-meanCandSMAP).^2/sum_xi_minus_xbar_squared));
%from Davis Statistics textbook pg 202

StandardDeviation = StandardError / (sqrt(1/length_Cotton_SheldonSz));

StandardError_newobs = sqrt(MSE*(1 + 1/length_Cotton_SheldonSz + (all_MAP-meanCandSMAP).^2/sum_xi_minus_xbar_squared));
%the variance of new indivual observation will be larger than the variance
%of means of observations (the variance of the means is related to the 
%standard error above). The variance in a new observation is what we are
%interested in


yhat_all = 5.673*all_MAP-269.87;

%Below plots data and regression line and confidence intervals
% figure(6)
% plot (Cotton_SheldonMAP, Cotton_SheldonSz, 'ko')
% hold on
% plot (all_MAP, yhat_all)
% % plot (all_MAP, yhat_all + StandardError, 'b--')
% % plot (all_MAP, yhat_all - StandardError, 'b--')
% % plot (all_MAP, yhat_all + StandardDeviation, 'b--')
% % plot (all_MAP, yhat_all - StandardDeviation, 'b--')
% plot (all_MAP, yhat_all + StandardError_newobs, 'r--')
% plot (all_MAP, yhat_all - StandardError_newobs, 'r--')
% xlabel('MAP (mm)')
% ylabel('S(z) (ppmV)')



%below plots the relative error versus absolute value of S(z)
%figure(7)
%plot (yhat_all, StandardError_newobs./yhat_all)
