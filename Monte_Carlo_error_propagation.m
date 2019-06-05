%This is a subroutine within PBUQ




%A warning is given the first time negative MAP or atmospheric CO2
%concentrations are calculated. Initializing these values now as zero and then changing
%their value if negative MAP or atm CO2 are calculated prevents multiple
%warnings from being printed.
warning_MAPzero = 0;
warning_MAPhigh = 0;
warningCa = 0;


numberofrandomvalues = 10000;
%The loop below selects random values for the variables that go into the
%paleobarometer from normal distributions with assigned standard 
%deviations 

for i = 1:numberofrandomvalues;
    
 
    %%%%%%%%%%%%%%%%%%%%%_____ d13Cs_____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    d13Cpc(:,i) = random('normal',bestd13Cpc, SD_d13Cc); %generates a normal pdf for d13Cpc
    temperature(:,i) = random('normal', MAT, SD_MAT); 
    %MAT is sometimes calculated from D47 T for input here (tempcalcmethod = 1).
    %in this case SD_MAT is actually SD_D47 T and the error associated
    %with the transfer function (MAT - carbonate formation T, i.e. the 
    %variable 'SE_MAATtocarbformTtransferfunctionerror') is set to equal zero.
    
    error_on_MAAT_to_Tcarbonateformation_transferfunction(:,i) = random('normal', 0, SE_MAATtocarbformTtransferfunctionerror);
    temperature(:,i) = 0.506 * temperature(:,i) + 17.974 + error_on_MAAT_to_Tcarbonateformation_transferfunction(:,i);
    
    d13Cs(:,i) = (d13Cpc(:,i)+1000)./((11.98-0.12*temperature(:,i))./1000+1)-1000;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%_____ S(z) _____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Szmethod == 1  % S(z) from MAP (Cotton and Sheldon 2012)
        
        if MAPmethod == 1 %MAP from depth to Bk
       
            %****X from Y error calculations********
            
            %random values are chosen from the 't' distribution for use in
            %the equations below to calculate condfidence intervals
            tstatistic = abs(random('t', length(log_d_to_Bk)-2));       

            standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
            %The variables D and H are used to calculate confidence intervals
            %on X values determined from measured Y values (Sokal and Rohlf
            %1981)
            D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
            %one Y for each X (as opposed to multiple Y values for each X in
            %the calibration)
        
            H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(log_d_to_Bk))...
                + ((log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk)).^2)/sum_of_squares_of_X));

            
            %The RHS of the equations below are the lower and upper confidence
            %limits for X (which is in this case MAP). They can be used to 
            %generate a probability density function for MAP. To make MAP 
            %array same length as other arrays, the lower confidence limits 
            %are calculated for the first half of the iterations and the 
            %upper confidence limits for the second half.
            if i < numberofrandomvalues/2
                mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D - H;   %the values of H and D change for each iteration because they are functions of tstatistic, for which a value is rondomly assinged from the t distribution each iteration
            else
                mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D + H;
            end
            
       
        elseif MAPmethod ==2 %MAP from CIA-K
            
            %****X from Y error calculations********
            
            %random values are chosen from the 't' distribution for use in
            %the equations below to calculate condfidence intervals
            tstatistic = abs(random('t', length(CIAK_Sheldon_calibration)-2));       

            standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
            %The variables D and H are used to calculate confidence intervals
            %on X values dtermined from measured Y values (Sokal and Rohlf
            %1981)
            D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
            %one Y for each X (as opposed to multiple Y values for each X in
            %the calibration)
        
            H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(CIAK_Sheldon_calibration))...
                + ((bestCIAK-mean(CIAK_Sheldon_calibration)).^2)/sum_of_squares_of_X));

            %see lines 68-72 for explanation of the code below
            if i < numberofrandomvalues/2
                ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D - H;
                mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
            else
                ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D + H;
                mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
            end
            
        elseif MAPmethod ==3
            mean_annual_precip(:,i) = random('normal', best_mean_annual_precip, SD_MAP);
        end 
         
        %This loop makes any negative or zero MAP values equal to 1 and
        %outputs a warning about invalid MAP values
        for j = 1:m-1  %m is the number of rows in the input matrix (number of different ages being considered). The first row consists of column headings
            if mean_annual_precip(j,i) <= 0
                mean_annual_precip(j,i) = 1;  %negative MAP values are not possible and result in imaginary values for atmospheric CO2 if MAP is used to calculate d13Cr
                if warning_MAPzero == 0;
                    sprintf('WARNING: negative MAP values calculated')
                    warning_MAPzero = 1;
                end    
            elseif mean_annual_precip(j,i) > 600
                if warning_MAPhigh == 0;
                    sprintf('WARNING: MAP values above calibrated range for S(z)')
                    warning_MAPhigh = 1;
                end 
            end
        end

       
        error_on_MAP_to_Sz_transferfunction (:,i) = random('normal', 0, SE_Sztransferfunctionerror);
    
        Sz(:,i) = (5.67*mean_annual_precip(:,i)-269.9)+ error_on_MAP_to_Sz_transferfunction (:,i);  %from Cotton and Sheldon 2012
    
    
    elseif Szmethod == 2 %S(z) from depth to Bk (Retallack 2009)
        error_on_depth_to_Bk_to_Sz_transferfunction(:,i) = random('normal', 0, SE_Sz_from_Bk_transferfunctionerror);
        %Sz(:,i) = 67.0*best_decompacted_depth_to_Bk+567 + error_on_depth_to_Bk_to_Sz_transferfunction(:,i);
        
        
        %revised to reflect minimum warm season S(z)
        Sz(:,i) = 35.3*best_decompacted_depth_to_Bk+588 + error_on_depth_to_Bk_to_Sz_transferfunction(:,i);

%Original PBUQ code for S(z) error below        
%         %****X from Y error calculations********
%         
%         %random values are chosen from the 't' distribution for use in
%         %the equations below to calculate condfidence intervals
%         tstatistic = abs(random('t', length(depth_to_Bk_Retallack_2009)-2));       
% 
%         standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
% 
%         %The variables D and H are used to calculate confidence intervals
%         %on X values dtermined from measured Y values (Sokal and Rohlf
%         %1981)
%         D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
% 
%         %one Y for each X (as opposed to multiple Y values for each X in
%         %the calibration)
% 
% 
%         H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(depth_to_Bk_Retallack_2009))...
%             + ((best_decompacted_depth_to_Bk-mean(depth_to_Bk_Retallack_2009)).^2)/sum_of_squares_of_X));
% 
%         %see lines 68-72 for explanation of the code below
%         if i < numberofrandomvalues/2
%             Sz(:,i) = mean(Sz_Retallack_2009) + b_yx * (best_decompacted_depth_to_Bk-mean(depth_to_Bk_Retallack_2009))/D - H;
%         else
%             Sz(:,i) = mean(Sz_Retallack_2009) + b_yx * (best_decompacted_depth_to_Bk-mean(depth_to_Bk_Retallack_2009))/D + H;
%         end
%end of original PBUQ code removed   
    
    elseif Szmethod == 3 %Montanez 2013 GCA soil order based S(z)
        
        for j = 1:m-1
            if soilorder(j) == 1
                Sz(j,i) = Mollisol_Sz(round(random('uniform',1,length(Mollisol_Sz))));
            elseif soilorder(j) == 2
                Sz(j,i) = Alfisol_Sz(round(random('uniform',1,length(Alfisol_Sz))));
            elseif soilorder(j) == 3
                Sz(j,i) = Aridisol_Sz(round(random('uniform',1,length(Aridisol_Sz))));
            elseif soilorder(j) == 4
                Sz(j,i) = Vertisol_Sz(round(random('uniform',1,length(Vertisol_Sz))));
            elseif soilorder(j) == 5
                Sz(j,i) = Andisol_Sz(round(random('uniform',1,length(Andisol_Sz))));
            elseif soilorder(j) == 6
                Sz(j,i) = Inceptisol_Sz(round(random('uniform',1,length(Inceptisol_Sz))));
            end
        end
        
    elseif Szmethod == 4 %other estimate of S(z)
         Sz(:,i) = random('normal', bestSz, SD_Sz);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%_____ d13Cr _____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if d13Crmethod ==1  %bulk paleosol OM
        d13Com(:,i) = random('normal',bestd13Com, SD_d13Com);
        OM_BigDelta_soil_paleosol(:,i) =random ('normal',best_OM_BigDelta_soil_paleosol,SD_OM_BigDelta_soil_paleosol);   
        %accounts for difference between paleosol bulk OM and soil OM; 
        %estimated for now, the magnitude of Big delta and its associated 
        %error should be quantified in the future.
        
        for j = 1:m-1
            bigdelta_respiredCO2_som(j,i) =random ('normal',best_bigdelta_respiredCO2_som(j),SD_d13Csom_to_d13Cr);   
            %accounts for difference between soil OM and soil-respired CO2; 
            %estimated for now, this error should be better quantified in the future
            d13Cr(j,i) = d13Com(j,i) + OM_BigDelta_soil_paleosol(i) + bigdelta_respiredCO2_som(j,i);
        end
        
    elseif d13Crmethod == 2    %occluded OM, the assumption here is that 
                               %there is no difference between d13C values of 
                               %occuluded OM and the OM in the soil before 
                               %diagenesis (i.e. carbonate perfectly protects OM)
        d13Com(:,i) = random('normal',bestd13Com, SD_d13Com);
        
        for j = 1:m-1
            bigdelta_respiredCO2_som(j,i) =random ('normal',best_bigdelta_respiredCO2_som(j),SD_d13Csom_to_d13Cr);   
            %estimated for now, this error should be quantified n the future
            d13Cr(j,i) = d13Com(j,i) + bigdelta_respiredCO2_som(j,i);
        end
    
    elseif d13Crmethod == 3    %d13Cr from d13Ca
        
        if d13Camethod == 1 %Tipple et al 2010
            d13Ca(:,i) = random('normal', bestd13Ca, SD_d13Ca);
        end
        
            
        if d13Camethod == 2 %Planktonic forams following Passey
            d13C_PF_carb(:,i) = random('normal', best_d13C_PF_carb, SD_d13C_PF_carb);
            epsilon13C_PF(:,i) = random('normal', best_epsilon13C_PF, 1.1); %Passey et al 2002 report one sigma = 1.1 for their estimate of epsilon
            d13Ca(:,i) = (d13C_PF_carb(:,i) + 1000)./(epsilon13C_PF(:,i)./1000 + 1) - 1000;
        end
        
        
        
        if d13Camethod ==3 %humid climate OM
        
        %****X from Y error calculations********
        
        d13Cplant(:,i) = random('normal', bestd13Cplant, SD_d13Cplant);    
        %This is not X from Y sections of code above because the equivalent 
        %variables (CIA-K or depth to Bk) is assumed to have no error
        
        
        %random values of are chosen from the 't' distribution for use in
        %the equations below to calculate confidence intervals
        tstatistic = abs(random('t', length(delta13Cplant_Arens_calibration)-2));       

        standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
        %The variables D and H are used to calculate confidence intervals
        %on X values dtermined from measured Y values (Sokal and Rohlf
        %1981)
        D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
        %one Y for each X (as opposed to multiple Y values for each X in
        %the calibration)

        %in order to calculate confidence interval around mean estimated X
        %value (equation for this is not given in Sokal and Rohlf) I removed the "1" added
        %to 1/n in the term multiplied by "D" just like Davis shows for
        %calculating standard errors of Y based on measurements of X (pg 204 and
        %202 of Davis texttbook)
        H = tstatistic/D * sqrt(s2YX*(D*(1/length(delta13Cplant_Arens_calibration))...
            + ((d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration)).^2)/sum_of_squares_of_X));

        %see lines 68-72 for explanation of the code below
           if i < numberofrandomvalues/2
               d13Ca(:,i) = mean(delta13Catm_Arens_calibration) + b_yx * (d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration))/D - H;
           else
               d13Ca(:,i) = mean(delta13Catm_Arens_calibration) + b_yx * (d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration))/D + H;
           end
        end   
        
        if d13Camethod == 4 %other d13Ca
            d13Ca(:,i) = random('normal', bestd13Ca, SD_d13Ca);
        end
        
        %Now we have d13Ca and want to calculate d13Cr
        %There is uncertainty associated with the magnitude of the
        %difference between d13C values of atmosphere CO2 and plant leaves
        %(d13Cp). There are two components of this uncertainty
        
        %One component of this uncertainty arises from uncertainty
        %in the value of MAP. The probability density function of MAP values
        %can be determined as follows:
        
        if MAP_or_no_MAP == 1
            
            if Szmethod > 1
                if MAPmethod == 1
                    
                    %****X from Y error calculations********
        
                    %random values are chosen from the 't' distribution for use in
                    %the equations below to calculate condfidence intervals
                    tstatistic = abs(random('t', length(log_d_to_Bk)-2));       

                    standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
                    %The variables D and H are used to calculate confidence intervals
                    %on X values determined from measured Y values (Sokal and Rohlf
                    %1981)
                    D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
                    %one Y for each X (as opposed to multiple Y values for each X in
                    %the calibration)

        
                    H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(log_d_to_Bk))...
                        + ((log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk)).^2)/sum_of_squares_of_X));

                    %see lines 68-72 for explanation of the code below
                    if i < numberofrandomvalues/2
                        mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D - H;
                    else
                        mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D + H;
                    end
                    
                elseif MAPmethod ==2 %MAP from CIA-K
                    
                    %****X from Y error calculations********
                    
                    %random values are chosen from the 't' distribution for use in
                    %the equations below to calculate condfidence intervals
                    tstatistic = abs(random('t', length(CIAK_Sheldon_calibration)-2));       

                    standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
                    %The variables D and H are used to calculate confidence intervals
                    %on X values dtermined from measured Y values (Sokal and Rohlf
                    %1981)
                    D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
                    %one Y for each X (as opposed to multiple Y values for each X in
                    %the calibration)

        
                    H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(CIAK_Sheldon_calibration))...
                        + ((bestCIAK-mean(CIAK_Sheldon_calibration)).^2)/sum_of_squares_of_X));

                    %see lines 68-72 for explanation of the code below
                    if i < numberofrandomvalues/2
                        ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D - H;
                        mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
                    else
                        ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D + H;
                        mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
                    end
                    
                elseif MAPmethod ==3
                    mean_annual_precip(:,i) = random('normal', best_mean_annual_precip, SD_MAP);
            
                end 
                
                %This loop makes any negative or zero MAP values equal to 1 and
                %outputs a warning about invalid MAP values
                for j = 1:m-1  %m is the number of rows in the input matrix (number of different ages being considered). The first row is the consists of column headings
                    if mean_annual_precip(j,i) <= 0
                        mean_annual_precip(j,i) = 1;  %negative MAP values are not possible and result in imaginary values for atmospheric CO2 if MAP is used to calculate d13Cr
                        if warning_MAPzero == 0;
                            sprintf('WARNING: negative MAP values calculated')
                            warning_MAPzero = 1;
                        end
                    end
                end
            end
        
            %And the probability density function of the magnitude of the carbon
            %isotope fractionation between leaf and atmosphere can be
            %determined by using the pdf of MAP values in the Kohn 2010 equation
            big_delta_leaf_atm(:,i) = 2.01 - 1.98 *10^-4 * altitude + 5.88 * log10(mean_annual_precip(:,i) + 300) + 0.00192 * abs(latitude); %Kohn 2010
        
        
        
            %A second component of the uncertainty associated with the magnitude of the
            %difference between d13C values of atmosphere CO2 and plant leaves
            %arises from uncertainty in Kohn's calibration. The uncertainty
            %associated with the calibration itself is calculated in the
            %subroutine d13Cr_error_calculator
            d13Cr_error_calculator 
            SE_BIG_DELTA = sqrt(MSE*(1/length_BIG_DELTA + (log10(best_mean_annual_precip+300)-mean_log_MAP_plus_300).^2/sum_xi_minus_xbar_squared));
            %Error for regression type 'Y from X.'
            %mean_log_MAP_plus_300 is related to the mean of MAP values over which the
            %calibration (Kohn 2010) has been made, whereas best_mean_annual_precip is
            %the best estimate for MAP for the paleosol of interest. The
            %standard error of the mean is used here, because we want to use mean
            %d13C values of plants in a whole ecosystem as an estimate of d13Cr, 
            %and we are not interested in the SE of individual new
            %observations (i.e. d13C values individual plants at a given
            %MAP)

            %The second component of the uncertainty can be considered as a normally 
            %distributed variable with a mean of zero. Randomly selected values 
            %from this distribution can be added to the values randomly selected
            %from the distribution of big deltas resulting from uncertain MAP values. This
            %is equivalent to randomly selecting values from a specified X-Y 
            %region on Kohn's plot of big delta versus MAP (i.e. first the MAP
            %range is selected and then the big delta range).
            Big_Delta_atm_leaf_error (:,i) = random('normal', 0,SE_BIG_DELTA);
            %This is using Kohn's calibration which relates Big Delta to MAP        
        
            %The resulting probability distribution of d13Cp can be determined
            %as follows:
            d13Cp(:,i) = (d13Ca(:,i) ...                         %incorporates error from uncertainty associated with d13Ca
                -(big_delta_leaf_atm(:,i)+...                      %incorporates error from uncertainty associated with MAP
                Big_Delta_atm_leaf_error(:,i)))./...              %incorporates error from Kohn (2010) calibration
                ((big_delta_leaf_atm(:,i)+...                     %incorporates error from uncertainty associated with MAP
                Big_Delta_atm_leaf_error(:,i))...                 %incorporates error from Kohn (2010) calibration
                ./1000 + 1);
                %Using Kohn et al 2010 calibration of d13Cp from MAP
        
        else %if no MAP value available then can still estimate d13Cplant from d13Ca, but uncertainty is larger
            
            big_delta_leaf_atm(:,i) = random('normal', best_big_delta_leaf_atm, 1);  %The value of 1 comes from the std deviation of means of 100 mm MAP bins (from 0-1000 mm) in Kohn 2010 compilation
            
            d13Cp(:,i) = (d13Ca(:,i)-big_delta_leaf_atm(:,i))./(big_delta_leaf_atm(:,i)./1000 + 1);
            %Using Kohn et al 2010 calibration of d13Cp from MAP
        end
        
        
        %Now we have d13Cplant and need to calculate d13Cr
        
        %There is uncertainty associated with the relationship (Bowling et al 2008) between the d13C
        %value of leaves (d13Cp) and the d13C value of soil respired CO2    
        Big_Delta_leaf_srCO2 (:,i) = random('normal', best_big_delta_leaf_srCO2,SD_BIG_DELTA_leaf_srCO2);
        
        %A probability density fucntion of d13Cr values can be determined
        %as follows:
        d13Cr(:,i) = (d13Cp(:,i)- Big_Delta_leaf_srCO2 (:,i))./(Big_Delta_leaf_srCO2 (:,i)/1000 + 1);   %Bowling et al 2008
              
   
    elseif d13Crmethod == 4   %other estimate of d13Cr
        d13Cr(:,i) = random('normal', bestd13Cr, SD_d13Cr);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%_____ d13Ca _____%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if d13Camethod == 1   %Tipple
        d13Ca(:,i) = random('normal', bestd13Ca, SD_d13Ca);
        
    elseif d13Camethod == 2  %planktonic foram carbonate
        d13C_PF_carb(:,i) = random('normal', best_d13C_PF_carb, SD_d13C_PF_carb);
        epsilon13C_PF(:,i) = random('normal', best_epsilon13C_PF, 1.1); %Passey et al 2002 report one sigma = 1.1 for their estimate of epsilon
        d13Ca(:,i) = (d13C_PF_carb(:,i) + 1000)./(epsilon13C_PF(:,i)./1000 + 1) - 1000;
        
    elseif d13Camethod == 3  %from d13C of well preserved OM using Arens et al 2000
        
        %****X from Y error calculations********
        
        d13Cplant(:,i) = random('normal', bestd13Cplant, SD_d13Cplant);
        
        %random values of are chosen from the 't' distribution for use in
        %the equations below to calculate condfidence intervals
        tstatistic = abs(random('t', length(delta13Cplant_Arens_calibration)-2));       

        standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485
        
        %The variables D and H are used to calculate confidence intervals
        %on X values dtermined from measured Y values (Sokal and Rohlf
        %1981)
        D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;
 
        %one Y for each X (as opposed to multiple Y values for each X in
        %the calibration)

        %in order to calculate confidence interval around mean estimated X value,
        %this equation is not given in Sokal and Rohlf, but I removed the "1" added
        %to 1/n in the term multipl;ied by "D" just like Davis shows for
        %calculating standard errors of Y based on measurements of X (pg 204 and
        %202 of Davis texttbook)
        H = tstatistic/D * sqrt(s2YX*(D*(1/length(delta13Cplant_Arens_calibration))...
            + ((d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration)).^2)/sum_of_squares_of_X));

        %see lines 68-72 for explanation of the code below
        if i < numberofrandomvalues/2
            d13Ca(:,i) = mean(delta13Catm_Arens_calibration) + b_yx * (d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration))/D - H;
        else
            d13Ca(:,i) = mean(delta13Catm_Arens_calibration) + b_yx * (d13Cplant(:,i)-mean(delta13Cplant_Arens_calibration))/D + H;
        end

        
        
    elseif d13Camethod == 4  %from d13Cr and MAP
        
        %There is uncertainty associated with relationship between the d13C
        %value of leaves (d13Cp) and the d13C value of soil respired CO2
        Big_Delta_leaf_srCO2 (:,i) = random('normal', best_big_delta_leaf_srCO2,SD_BIG_DELTA_leaf_srCO2); %Bowling et al 2008
        
        
        
        %Therefore d13Cp has some degree of uncertainty. The probability distribution
        %of d13Cp values can de determined as follows:
        d13Cp(:,i) = d13Cr(:,i).*(Big_Delta_leaf_srCO2 (:,i)./1000+1) + Big_Delta_leaf_srCO2 (:,i);  
        %Uses Bowling et al 2008 to determine d13Cp from d13Cr
        
        
        
        
        %There is also uncertainty associated with the magnitude of the
        %difference between d13C values of atmosphere CO2 and plant leaves
        %(d13Cp). There are two components of this uncertainty
        
        %One component of this uncertainty arises from uncertainty
        %in the value of MAP. The probability distribution of MAP values
        %can be determined as follows:
        if Szmethod > 1  %If Szmethod equals 2,3 or 4, then a pdf for MAP has not yet been defined
            if MAPmethod == 1
                
                %****X from Y error calculations********

                %random values are chosen from the 't' distribution for use in
                %the equations below to calculate condfidence intervals
                tstatistic = abs(random('t', length(log_d_to_Bk)-2));       

                standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485

                %The variables D and H are used to calculate confidence intervals
                %on X values dtermined from measured Y values (Sokal and Rohlf
                %1981)
                D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;

                %one Y for each X (as opposed to multiple Y values for each X in
                %the calibration)

                H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(log_d_to_Bk))...
                    + ((log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk)).^2)/sum_of_squares_of_X));

                %see lines 68-72 for explanation of the code below
                if i < numberofrandomvalues/2
                    mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D - H;
                else
                    mean_annual_precip(:,i) = mean(MAP_R_2005) + b_yx * (log_best_decompacted_depth_to_Bk-mean(log_d_to_Bk))/D + H;
                end

            elseif MAPmethod ==2 %MAP from CIA-K
                
                %****X from Y error calculations********
                
                %random values are chosen from the 't' distribution for use in
                %the equations below to calculate condfidence intervals
                tstatistic = abs(random('t', length(CIAK_Sheldon_calibration)-2));       

                standard_error_regression_coefficient = sqrt(s2YX/sum_of_squares_of_X); %Sokal & Rohlf, pg 485

                %The variables D and H are used to calculate confidence intervals
                %on X values dtermined from measured Y values (Sokal and Rohlf
                %1981)
                D = b_yx^2 -  tstatistic^2 * standard_error_regression_coefficient^2;

                %one Y for each X (as opposed to multiple Y values for each X in
                %the calibration)

                H = tstatistic/D * sqrt(s2YX*(D*(1+1/length(CIAK_Sheldon_calibration))...
                    + ((bestCIAK-mean(CIAK_Sheldon_calibration)).^2)/sum_of_squares_of_X));

                %see lines 68-72 for explanation of the code below
                if i < numberofrandomvalues/2
                    ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D - H;
                    mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
                else
                    ln_mean_annual_precip(:,i) = mean(ln_MAP_Sheldon_calibration) + b_yx * (bestCIAK-mean(CIAK_Sheldon_calibration))/D + H;
                    mean_annual_precip(:,i) = exp(ln_mean_annual_precip(:,i));
                end

            elseif MAPmethod ==3
                mean_annual_precip(:,i) = random('normal', best_mean_annual_precip, SD_MAP);

            end 

            %This loop makes any negative or zero MAP values equal to 1 and
            %outputs a warning about invalid MAP values
            for j = 1:m-1  %m is the number of rows in the input matrix (number of different ages being considered). The first row is the consists of column headings
                if mean_annual_precip(j,i) <= 0
                    mean_annual_precip(j,i) = 1;  %negative MAP values are not possible and result in imaginary values for atmospheric CO2 if MAP is used to calculate d13Cr
                    if warning_MAPzero == 0;
                        sprintf('WARNING: negative MAP values calculated')
                        warning_MAPzero = 1;
                    end
                end
            end
        end
        
        %And the probability distribution of the magnitude of the carbon
        %isotope fractionation between leaf and atmosphere can be
        %determined by using the pdf of MAP values in Kohn 2010 equation
        big_delta_leaf_atm(:,i) = 2.01 - 1.98 *10^-4 * altitude + 5.88 * log10(mean_annual_precip(:,i) + 300) + 0.00192 * abs(latitude); %Kohn 2010
        
        
        
        %A second component of the uncertainty associated with the magnitude of the
        %difference between d13C values of atmosphere CO2 and plant leaves
        %arises from uncertainty in Kohn's calibration. The uncertainty
        %associated with the calibration itself is calculated in the
        %subroutine d13Cr_error_calculator
        d13Cr_error_calculator 
        SE_BIG_DELTA = sqrt(MSE*(1/length_BIG_DELTA + (log10(best_mean_annual_precip+300)-mean_log_MAP_plus_300).^2/sum_xi_minus_xbar_squared));
            %Error for regression type 'Y from X.'
            %mean_log_MAP_plus_300 is related to the mean of MAP values over which the
            %calibration (Kohn 2010) has been made, whereas best_mean_annual_precip is
            %the best estimate for MAP for the paleosol of interest. The
            %standard error of the mean is used here, because we want to use mean
            %d13C values of plants in a whole ecosystem as an estimate of d13Cr, 
            %and we are not interested in the SE of individual new
            %observations (i.e. d13C values individual plants at a given
            %MAP)
        
        %The second component of the uncertainty can be as a normally 
        %distributed variable with a mean of zero. Randomly selected values 
        %from this distribution can be added to the values randomly selected
        %from the distribution of resulting from uncertain MAP values. This
        %is equivalent to randomly selecting values from a specified X-Y 
        %region on Kohn's plot of big delta versus MAP (i.e. first the MAP
        %range is selected and then the big delta range).
        Big_Delta_atm_leaf_error (:,i) = random('normal', 0,SE_BIG_DELTA);
        %This is using Kohn's calibration which relates Big Delta to MAP
        
        
        %The resulting probability density function of d13Ca can be determined
        %as follows:
        d13Ca(:,i) = d13Cp(:,i).*...                         %incorporates error from Bowling et al 2008 relationship and from uncertainty associated with d13Cr
            ((big_delta_leaf_atm(:,i) + ...                  %incorporates error from uncertainty associated with MAP
            Big_Delta_atm_leaf_error(:,i))...                %incorporates error from Kohn (2010) calibration
            ./1000+1) + ...    
            (big_delta_leaf_atm(:,i) + ...                   %incorporates error from uncertainty associated with MAP
            Big_Delta_atm_leaf_error(:,i));                  %incorporates error from Kohn (2010) calibration
             %Uses Kohn 2010 to determine d13Ca from d13Cp
             
             
    
             
    elseif d13Camethod == 5   %other estimate of d13Ca
        d13Ca(:,i) = random('normal', bestd13Ca, SD_d13Ca);
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%end of d13Ca section%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    
%%%%%%%%%%%%%___EQUATION FOR PALEOSOL BAROMETER____%%%%%%%%%%%%%%%%%%%%%%%%

%this generates a pdf of atmospheric CO2
%concentration that estimates overall error

    SminusR (:,i) = (d13Cs(:,i)-1.0044*d13Cr(:,i)-4.4);   %numerator of barometer equation
    AminusS (:,i) = (d13Ca(:,i)-d13Cs(:,i));              %denominator of barometer equation
    SminusR_divby_AminusS (:,i) = (d13Cs(:,i)-1.0044*d13Cr(:,i)-4.4)./(d13Ca(:,i)-d13Cs(:,i));      
    %slope of equation (in atm CO2 versus S(z) space), this is 'R'
    
    Ca(:,i) = Sz(:,i).*(d13Cs(:,i)-1.0044*d13Cr(:,i)-4.4)./(d13Ca(:,i)-d13Cs(:,i));                     
    %atmopsheric CO2 concentration
    
    for j = 1:m-1  %m is the number of rows in the input matrix (number of different ages being considered). The first row is the consists of column headings
        if Ca(j,i) <= 0
            Ca(j,i) = 1;  %negative Ca values are not possible. If Ca < 0, make Ca = 1
            if warningCa == 0;
                sprintf('WARNING: negative Ca values calculated')
                warningCa = 1;
            end
        end
    end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end


Ca=Ca';
SminusR_divby_AminusS = SminusR_divby_AminusS';


atm_CO2_estimate (1,:) = median(Ca)   %median of a matrix gives mean of each column
atm_CO2_estimate (2,:) = std(Ca);
atm_CO2_estimate= atm_CO2_estimate';   %puts median for each age in first column and 1 sigma in second column
max_atm_CO2 = mean(Ca) + std (Ca);
max_atm_CO2 = max_atm_CO2';
RATIO = mean(SminusR_divby_AminusS);        
RATIO = RATIO';


figure(1)  % plot histograms of delta values
hist(d13Cs(1,:),50)    %change the index here to plot histograms for paleosols in rows other than 1
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','r')
hold on
hist(d13Ca(1,:),50)  %change the index here to plot histograms for paleosols in rows other than 1
hist(1.0044*d13Cr(1,:)+4.4, 50)  %change the index here to plot histograms for paleosols in rows other than 1

figure(2) % plot histogram of S(z) values
hist(Sz(1,:),50)   %change the index here to plot histograms for paleosols in rows other than 1


% Tanomaly = 4;
% numberofpaleosolsaveraged = [1:1:50];
% length_numberofpaleosolsaveraged = length(numberofpaleosolsaveraged);
% 
% for j =1:length_numberofpaleosolsaveraged
%     
%     for i=1:1000
%         %randomly select values from distributions of Ca
%         X = round(random('Uniform', 1,length(Ca),numberofpaleosolsaveraged(j),1));  %randomly select index of Ca matrix (number of random index numbers is, for each iteration, is equal to the number of paleosols averaged)
%         for k = 1:j 
%             CO2_high(k) = Ca(X(k),12);     %12 is the SJII
%             CO2_low(k) = Ca(X(k),11);      %11 is Nebraskan Mollisol
%         end
%         %Now average the randomly selected Ca values- this approximates averaging different paleosols (assuming that secular variation in CO2 is insignificant)
%         CO2_high_median(i,j) = median(CO2_high);    %such means are calculated 'i' number of times to account for the uncertainty of individual estimates (i.e. a monte carlo simulations of mean Ca for a high CO2 and low Co2 estimates)
%         CO2_low_median(i,j) = median(CO2_low);
%         ESS(i,j) = Tanomaly/(log2(CO2_high_median(i,j)/CO2_low_median(i,j)));
%     end
%    
% end
% 
% medianESS = median(ESS);
% sortedESS = sort(ESS);
% upper95percent_ESS=sortedESS(round(0.95*length(sortedESS)),:);  %95th percntile value
% lower5percent_ESS=sortedESS(round(0.05*length(sortedESS)),:);  %5th percentile values
% 
% figure(10)
% errorbar(numberofpaleosolsaveraged, medianESS, medianESS-lower5percent_ESS, upper95percent_ESS-medianESS, 'o-r');
% xlabel('number of paleosols averaged')
% ylabel('ESS')