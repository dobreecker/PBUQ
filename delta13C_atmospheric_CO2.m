%This is a subroutine within PBUQ

%The carbon isotope composition of atmospheric CO2 can be determined in a
%number of ways: it can be calcuated from the d13C values of contemporaneous
%organic matter or marine carbonates. The d13C value of atmospheric CO2 is 
%closely related to the d13C value of respired CO2. These two values are expected to
%covary and in practice, one is typically calculated from another. PBUQ 
%incorporates the uncertainty associated with this covariation
%and provides numerous options for determining the d13C value of
%atmospheric CO2.


%If d13Camethod = 0, then d13Ca has not yet been assigned and therefore
%needs to be assigned here
if d13Camethod == 0
    
    %ask user to select the method preferred for d13Ca
    d13Camethod = input('How do you want to determine the d13C value of atmospheric CO2? \n 1) using Tipple et al 2010 \n 2) contemporaneous planktonic foram carbonates (following Passey et al 2002) \n 3) from the d13C value of well-preserved, humid climate OM (using Arens et al 2000) \n 4) from MAP and d13Cr (using Kohn 2010) \n 5) other d13Ca estimate \n')
        
    if d13Camethod == 1 %using Tipple
        
        
        load tipple_et_al_2010.txt
        
        for k = 1:m-1
            for j= 1:length(tipple_et_al_2010)
                if age(k) < tipple_et_al_2010(j,1)   %find ages in Tipple that correspond to ages of paleosols. Ages are in column 1 of tipple_et_al_2010.txt
                    bestd13Ca(k) = mean(tipple_et_al_2010(j-2:j+1,2))';    %mean of five local d13Ca values
                    stdevd13Ca(k) = std(tipple_et_al_2010(j-2:j+1,2))';    %standard deviation of five local d13Ca values (error is not propagated through the calculation of individual d13Ca estimates)
                break
                end
            end
        end
        
        bestd13Ca=bestd13Ca';
        SD_d13Ca=stdevd13Ca';
        
    elseif d13Camethod == 2 %contemp. marine carb (Passey et al 2002)
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 21 %finds the column for contemporaneous marine planktonic foram carbonate d13C values
                best_d13C_PF_carb = paleosol_CO2_error_quant(2:m, i);
            end
        end
        best_epsilon13C_PF = 7.9;                           %from Passey et al 2002
        bestd13Ca = (best_d13C_PF_carb + 1000)./(best_epsilon13C_PF./1000 + 1) - 1000;
        
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 23 %finds the column for error associated with measured d13C values of marine planktonic foram carbonate
                SD_d13C_PF_carb = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        
    elseif d13Camethod == 3 %from well preserved, humid climate OM (Arens et al 2000)
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 17  %finds the column for undecomposed bulk organic matter
                bestd13Cplant = paleosol_CO2_error_quant(2:m, i);
                bestd13Ca = (paleosol_CO2_error_quant(2:m, i) + 18.67)./1.1; 
                 %This equation is from Arens et al 2000 and works if OM formed in a climate that
                %was humid enough that d13Com is not a function of MAP.
                %ALSO, THIS EQUATION IS NOT APPROPRIATE FOR GRASSES!
                %THIS APPROACH WORKS BEST FOR WHOLE ECOSYSTEM ORGANIC CARBON, NOT
                %INDIVIDUAL SPECIES (such as invidual roots!). The error
                %caluclated assumes the d13C value of om is an ecosystem mean.
            end
        end
        
        for i=1:n 
            if paleosol_CO2_error_quant(1,i) == 20
                SD_d13Cplant = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %***********X FROM Y************ 
        %(see lines 57-63 of subroutine soil_derived_component_of_soil_CO2 
        %for explanation)
        
        load Arensetal_2000_data.txt      %load calibration data

        delta13Catm_Arens_calibration = Arensetal_2000_data(:,1);
        delta13Cplant_Arens_calibration = Arensetal_2000_data(:,2);

        %Arens et al calculated 99% confidence intervals using method described in
        %Sokal and Rohlf. Their confodence intervals are small because they used
        %'ecosystem average' d13Cp values, which are means of a number of different
        %individuals and plant types. Therefore, confidence intervals are
        %around the mean rather than individual observations

        %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
        sum_of_squares_of_Y = sum((delta13Cplant_Arens_calibration-mean(delta13Cplant_Arens_calibration)).^2); %pg 55
        sum_of_squares_of_X = sum((delta13Catm_Arens_calibration-mean(delta13Catm_Arens_calibration)).^2);     %pg 55
        sum_xy = sum((delta13Catm_Arens_calibration-mean(delta13Catm_Arens_calibration)).*(delta13Cplant_Arens_calibration-mean(delta13Cplant_Arens_calibration)));  %Eq 14.4 page 467

        b_yx = sum_xy/sum_of_squares_of_X; %slope of regression line

        Y_intercept = mean(delta13Cplant_Arens_calibration)-b_yx*mean(delta13Catm_Arens_calibration);

        yhat = b_yx * delta13Catm_Arens_calibration + Y_intercept;

        unexplained_sum_of_squares = sum((delta13Cplant_Arens_calibration-yhat).^2);  %pg 469, also called error sum of squares

        s2YX = unexplained_sum_of_squares/(length(delta13Cplant_Arens_calibration)-2);   %pg 473, same as MSE above from Davis. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
      
    elseif d13Camethod == 4 %If user wants to calculate d13Ca using MAP, then need to know how user wants to determine MAP
        
        %beginning of MAP section of code
        if Szmethod > 1  %If Szmethod equals 2,3 or 4, then best_mean_annual_precip has not yet been defined
            
            %ask user to select the method preferred for MAP
            MAPmethod = input ('How do you want to determine MAP for calculating d13Ca? \n 1) from depth to B(k) (Retallack 2005) \n 2) from CIA-K (Sheldon et al 2002) \n 3) other estimate of MAP \n')
    
            if MAPmethod == 1 %MAP from depth to Bk
                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 8
                        best_depth_to_Bk = paleosol_CO2_error_quant(2:m, i);
                    end
                end

                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 9
                        best_burial_depth = paleosol_CO2_error_quant(2:m, i);
                    end
                end

                %error associated with depth to Bk and burial depth are not
                %considered here, but when propagated are likely smaller than other
                %errors
                best_decompacted_depth_to_Bk = best_depth_to_Bk./(-0.62./(0.38./exp(0.17*best_burial_depth)-1));
                %equation for decompaction for Aridisols (Sheldon and Tabor 2009)
                
                log_best_decompacted_depth_to_Bk = log(best_decompacted_depth_to_Bk);
                %log transformation used before fitting regression line
                
                best_mean_annual_precip = (log_best_decompacted_depth_to_Bk-2.97)/.0021; %This is the best fit regression line if MAP is considered the independent variable.  If depth to Bk considered the independent variable, 137.24+6.45*best_decompacted_depth_to_Bk-0.013*best_decompacted_depth_to_Bk.^2 is best fit, from Retallack 2005, who reported SE = 147 mm
                %This is the best fit regression line if MAP is considered the 
                %independent variable.  If depth to Bk considered the independent 
                %variable then the best fit polynomial is 
                %137.24 + 6.45*best_decompacted_depth_to_Bk - 
                %0.013*best_decompacted_depth_to_Bk.^2 (Retallack 2005)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %***********X FROM Y************ 
                %(see lines 57-63 of subroutine soil_derived_component_of_soil_CO2 
                %for explanation)
                
                load Retallack_2005_data.txt     %load calibration data
                

                MAP_R_2005 = Retallack_2005_data(:,1);
                depth_to_Bk_R_2005 = Retallack_2005_data(:,2);
                log_d_to_Bk = log(depth_to_Bk_R_2005); %natural log

                %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
                sum_of_squares_of_Y = sum((log_d_to_Bk-mean(log_d_to_Bk)).^2); %pg 55
                sum_of_squares_of_X = sum((MAP_R_2005-mean(MAP_R_2005)).^2);     %pg 55
                sum_xy = sum((MAP_R_2005-mean(MAP_R_2005)).*(log_d_to_Bk-mean(log_d_to_Bk)));  %Eq 14.4 page 467

                b_yx = sum_xy/sum_of_squares_of_X; %slope of regression line

                Y_intercept = mean(log_d_to_Bk)-b_yx*mean(MAP_R_2005);

                yhat = b_yx * MAP_R_2005 + Y_intercept;

                unexplained_sum_of_squares = sum((log_d_to_Bk-yhat).^2);  %pg 469, also called error sum of squares

                s2YX = unexplained_sum_of_squares/(length(log_d_to_Bk)-2);   %pg 473, same as MSE above from Davis. 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
            elseif MAPmethod == 2   %MAP from CIA-K
                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 10 % finds column with CIA-K values
                        bestCIAK = paleosol_CO2_error_quant(2:m, i);
                        best_mean_annual_precip = exp((bestCIAK+172.75)/35.7914); %This is the best fit line if MAP in regression where MAP is considered as the independent variable. If CIA-K is considered the independent variable, then the best fit line is: 221.1*exp(0.02*paleosol_CO2_error_quant(2:m, i)) Sheldon et al 2002
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %***********X FROM Y************ 
                    %(see lines 57-63 of subroutine soil_derived_component_of_soil_CO2 
                    %for explanation)
                    
                    load Sheldonetal_2002_CIAKdata.txt      %load calibration data

                    ln_MAP_Sheldon_calibration = Sheldonetal_2002_CIAKdata(:,1);
                    CIAK_Sheldon_calibration = Sheldonetal_2002_CIAKdata(:,2);


                    %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
                    sum_of_squares_of_Y = sum((CIAK_Sheldon_calibration-mean(CIAK_Sheldon_calibration)).^2); %pg 55
                    sum_of_squares_of_X = sum((ln_MAP_Sheldon_calibration-mean(ln_MAP_Sheldon_calibration)).^2);     %pg 55
                    sum_xy = sum((ln_MAP_Sheldon_calibration-mean(ln_MAP_Sheldon_calibration)).*(CIAK_Sheldon_calibration-mean(CIAK_Sheldon_calibration)));  %Eq 14.4 page 467

                    b_yx = sum_xy/sum_of_squares_of_X; %slope of regression line

                    Y_intercept = mean(CIAK_Sheldon_calibration)-b_yx*mean(ln_MAP_Sheldon_calibration);

                    yhat = b_yx * ln_MAP_Sheldon_calibration + Y_intercept;

                    unexplained_sum_of_squares = sum((CIAK_Sheldon_calibration-yhat).^2);  %pg 469, also called error sum of squares, I think, identical to SSE above from Davis

                    s2YX = unexplained_sum_of_squares/(length(CIAK_Sheldon_calibration)-2);   %pg 473, same as MSE above from Davis. The denominator if degrees of freedom, which may not be calculated correctly for the case of multiple Y values for each value of X- I may need to substract 2 from 'a', the number of groups in the anova.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 end
        
            elseif MAPmethod ==3   %other estimate of MAP
                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 11
                        best_mean_annual_precip = paleosol_CO2_error_quant(2:m, i);
                    end
                end
        
                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 12
                        SD_MAP = paleosol_CO2_error_quant(2:m, i);
                    end
                end
            end
        end
        %end of MAP section of code
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %altitude and latitude needed for Kohn's 2010 equation
        for i=1:n 
            if paleosol_CO2_error_quant(1,i) == 24
                latitude = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        for i=1:n 
            if paleosol_CO2_error_quant(1,i) == 25
                altitude = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        
        best_big_delta_leaf_atm = 2.01 - 1.98 *10^-4 * altitude + 5.88 * log10(best_mean_annual_precip + 300) + 0.00192 * abs(latitude); %Kohn 2010
        bestd13Cp = bestd13Cr.*(best_big_delta_leaf_srCO2./1000+1) + best_big_delta_leaf_srCO2;  
        %Uses Bowling et al to determine d13C of leaves from d13Cr
        
        bestd13Ca = bestd13Cp.*(best_big_delta_leaf_atm./1000+1) + best_big_delta_leaf_atm;     
        
        
        %d13Ca error is calculated in MC subroutine for d13Camethod = 4
        
    elseif d13Camethod == 5 %d13C values of atmospheric CO2 estimated some other way
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 22
                bestd13Ca = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 23  
                SD_d13Ca = paleosol_CO2_error_quant(2:m, i);
            end
        end
    end  
end