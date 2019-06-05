%This is a subroutine within PBUQ

%The carbon isotope composition of soil respired CO2 can be determined in a
%number of way: it can be calcuated from the d13C values of contemporaneous
%organic matter, atmospheric CO2, or marine carbonates (by way of atmospheric CO2). 
%The uncertainty associated with determinations of d13Cr is currently not 
%well quantified, largely because the d13C values of ecosystem carbon pools 
%and fluxes in calcic soils have not been sufficiently studied.


%Initialize some values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the d13C value of leaves is typically lower than the d13C value of soil
%respired CO2 (for forests)
best_big_delta_leaf_srCO2 = -2.6; %from Bowling et al 2008 
%and this difference has an associated uncertainty
SD_BIG_DELTA_leaf_srCO2 = 0.6; %estimated from figure 4 of Bowling et al 2008

%The  d13C value of bulk paleosol OM may not be the same as the d13C value
%of primary OM and therefore primary OM should be measured if possible
%(i.e. OM occluded in paleosol carbonate). In the absence of better 
%information, the difference between bulk paleosol and primary OM is set to 
%zero here, but can be adjusted in the future if determined appropriate.
best_OM_BigDelta_soil_paleosol = 0;  
%Below is an estimate of the error associated with the difference between 
%the d13C value of bulk paleosol OM and the d13C value of primary OM
SD_OM_BigDelta_soil_paleosol = 0.5; 

%The d13C value of soil repsired CO2 probably does not equal the d13C value
%of primary OM. 
best_d13Cr_minus_d13Com_A_horizon = 0.5; %based on values compiled by Bowling 2008
best_d13Cr_minus_d13Com_B_horizon = -1.0; 
%based on the value above and the difference between d13C values of OM in 
%A and B horizon of am archived soil (Torn et al 2002).

SD_d13Csom_to_d13Cr = 0.5;
%An estimate of the error

d13Camethod = 0; %Some options below involve calculating d13Ca first and then d13Cr.
%If d13Ca is calculated in this subroutine, then it does not need to be calculated in the next
%subroutine, and this value becomes nonzero, which allows the d13Ca subroutine to be skipped.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%ask user to select the method preferred for d13Cr
d13Crmethod = input ('How do you want to determine the d13C value of respired CO2? \n 1) bulk paleosol organic matter \n 2) OM occluded in paleosol carbonate or otherwise well preserved in paleosol \n 3) from d13Ca \n 4) other estimate of d13C value of respired CO2 \n')
 

if d13Crmethod == 1
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 15 %finds column for d13C values of bulk paleosol organic matter
            bestd13Com = paleosol_CO2_error_quant(2:m, i);
        end
    end
     
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 19 %finds column for standard deviations of d13C values of paleosol organic matter (in this case bulk)
            SD_d13Com = paleosol_CO2_error_quant(2:m, i);  
            %This is analytical error, error in determining d13Cr from 
            %d13Com is calculated in MC subroutine
        end
    end
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 27 %finds column with the numbers 1 or 2 designating A or B horizon, respectively
            om_AorB = paleosol_CO2_error_quant(2:m, i);  
        end
    end
    
    %the loop below assigns magnitudes for the difference between d13C
    %values of soil respired CO2 and om according to the horizon from which
    %the OM was collected and then calculates best d13Cr values.
    for j = 1:m-1    %there are m-1 rows with paleosol data in input text file
        if om_AorB(j) == 1
            best_bigdelta_respiredCO2_som(j,1) = best_d13Cr_minus_d13Com_A_horizon;
        else
            best_bigdelta_respiredCO2_som(j,1) = best_d13Cr_minus_d13Com_B_horizon;
        end
        
        bestd13Cr(j,1) = bestd13Com(j) + best_OM_BigDelta_soil_paleosol + best_bigdelta_respiredCO2_som(j,1);
    end
            
    
    %d13Cr error is calculated in MC subroutine for d13Crmethod = 1    
    
elseif d13Crmethod == 2
    for i=1:n
        if paleosol_CO2_error_quant(1,i) == 16 %finds column for d13C values of OM occluded in paleosol carbonate
            bestd13Com = paleosol_CO2_error_quant(2:m, i);
        end
    end
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 19 %finds column for standard deviations of d13C values of organic matter
            SD_d13Com = paleosol_CO2_error_quant(2:m, i); %This is analytical error, error in determining d13Cr from d13Com is calculated in MC subroutine
        end
    end
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 27 %finds column with 1s or 2s designating A or B horizon, respectively
            om_AorB = paleosol_CO2_error_quant(2:m, i);  
        end
    end
    
    %the loop below assigns magnitudes for the difference between d13C
    %values of soil respired CO2 and om according to the horizon from which
    %the OM was collected.
    for j = 1:m-1
        if om_AorB(j) == 1
            best_bigdelta_respiredCO2_som(j,1) = best_d13Cr_minus_d13Com_A_horizon;
        else
            best_bigdelta_respiredCO2_som(j,1) = best_d13Cr_minus_d13Com_B_horizon;
        end
        
        bestd13Cr(j,1) = bestd13Com(j) + best_bigdelta_respiredCO2_som(j,1);
    end
    %d13Cr error is calculated in MC subroutine for d13Crmethod = 2
    
    

     %method 4 before 3 because method 4 is short adn method 3 is quite long   
elseif d13Crmethod == 4
    for i=1:n
        if paleosol_CO2_error_quant(1,i) == 18
            bestd13Cr = paleosol_CO2_error_quant(2:m, i);
        end
    end

    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 19 %finds column for standard deviations of d13C values of respired CO2 estimated some other way
            SD_d13Cr = paleosol_CO2_error_quant(2:m, i);
        end
    end




elseif d13Crmethod == 3 
%If user wants to calculate d13Cr from d13Ca, then need to know how user 
%wants to determine d13Ca and MAP (if MAP not already calculated). 
%There is also an option if no MAP estimate is available, but uncertainty 
%is higher for this option    

    %First determine d13Ca%%%%%%%%%%%%%%%%%
    
    %ask user to select the method preferred for d13Ca
    d13Camethod = input('How do you want to determine the d13C value of atmospheric CO2? \n 1) using Tipple et al 2010 \n 2) contemporaneous phytoplankton carbonates (following Passey et al 2002) \n 3) from the d13C value of well-preserved, humid climate OM (using Arens et al 2000) \n 4) other d13Ca estimate \n')
        
    if d13Camethod == 1   %from Tipple et al (2010)
        
        
        load tipple_et_al_2010.txt
        
        
        for k = 1:m-1
            for j= 1:length(tipple_et_al_2010)
                if age(k) < tipple_et_al_2010(j,1)
                    %find ages in Tipple that correspond to ages of paleosols. 
                    %Ages are in column 1 of tipple_et_al_2010.txt  
                    
                    bestd13Ca(k) = mean(tipple_et_al_2010(j-2:j+1,2));     
                    %mean of five local d13Ca values
                    stdevd13Ca(k) = std(tipple_et_al_2010(j-2:j+1,2));     
                    %standard deviation of five local d13Ca values 
                    %(error is not propagated through the calculation of 
                    %individual d13Ca estimates) 
                    
                break
                end
            end
        end
        
        bestd13Ca=bestd13Ca';
        SD_d13Ca=stdevd13Ca';
                  
    elseif d13Camethod == 2
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 21  %finds column with d13C values of contemporaneous marine planktonic foram carbonate
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
    
     elseif d13Camethod == 3 %humid climate OM using Arens et al 2000
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 17  %finds the column for undecomposed bulk organic matter from penecomtempraneous deposit
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
        %Sokal and Rohlf. Their confidence intervals are small because they used
        %'ecosystem average' d13Cp values, which are means of a number of different
        %individuals and plant types. Therefore, confidence intervals are
        %around the mean rather than new observations.

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
        
        
    elseif d13Camethod == 4
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 22  %finds column with d13C values of atmospheric CO2 estimated some other way
                bestd13Ca = paleosol_CO2_error_quant(2:m, i);
            end
        end
        
        for i=1:n
            if paleosol_CO2_error_quant(1,i) == 23  %finds column with error associated with d13C values of atmospheric CO2 estimated some other way
                SD_d13Ca = paleosol_CO2_error_quant(2:m, i);
            end
        end
    end
    %end of d13Ca subsection%%%%%%%%%%%%
    
    
    
    
    
    %Now determine MAP for calculating d13Cr from d13Ca. First ask user if
    %MAP should be used
    MAP_or_no_MAP = input('Do you want to use an estimate of MAP to calculate d13Cr from d13Ca using Kohn 2010? \n 1) Yes \n 2) No \n')
    
    if MAP_or_no_MAP == 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %best values for MAP may or may not have already been assigned
        %If Szmethod equals 2, 3 or 4, then best_mean_annual_precip has not yet been 
        %defined and therefore needs to be
        if Szmethod > 1 
            
            %below is the same as the code in subroutine
            %'soil_derived_component_of_soil_CO2'
        
            MAPmethod = input ('How do you want to determine MAP for calculating d13Cr? \n 1) from depth to B(k) (Retallack 2005) \n 2) from CIA-K (Sheldon et al 2002) \n 3) other estimate of MAP \n')
    
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
    
                best_decompacted_depth_to_Bk = best_depth_to_Bk./(-0.62./(0.38./exp(0.17*best_burial_depth)-1));
                %equation for decompaction for Aridisols (Sheldon and Tabor 2009)
                
                log_best_decompacted_depth_to_Bk = log(best_decompacted_depth_to_Bk);
                %log transformation used before fitting regression line
                
                best_mean_annual_precip = (log_best_decompacted_depth_to_Bk-2.97)/.0021; %This is the best fit regression line if MAP is considered the independent variable.  If depth to Bk considered the independent variable, 137.24+6.45*best_decompacted_depth_to_Bk-0.013*best_decompacted_depth_to_Bk.^2 is best fit, from Retallack 2005
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
                log_d_to_Bk = log(depth_to_Bk_R_2005); %natural log (log transformation of depth tp Bk for regression)


                %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
                sum_of_squares_of_Y = sum((log_d_to_Bk-mean(log_d_to_Bk)).^2); %pg 55
                sum_of_squares_of_X = sum((MAP_R_2005-mean(MAP_R_2005)).^2);     %pg 55
                sum_xy = sum((MAP_R_2005-mean(MAP_R_2005)).*(log_d_to_Bk-mean(log_d_to_Bk)));  %Eq 14.4 page 467

                b_yx = sum_xy/sum_of_squares_of_X; %slope of regression line

                Y_intercept = mean(log_d_to_Bk)-b_yx*mean(MAP_R_2005);

                yhat = b_yx * MAP_R_2005 + Y_intercept;

                unexplained_sum_of_squares = sum((log_d_to_Bk-yhat).^2);  %pg 469, also called error sum of squares, 

                s2YX = unexplained_sum_of_squares/(length(log_d_to_Bk)-2);   %pg 473, same as MSE above from Davis. 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                
            elseif MAPmethod == 2   %MAP from CIA-K
                for i=1:n 
                    if paleosol_CO2_error_quant(1,i) == 10 % finds column with CIA-K values
                        bestCIAK = paleosol_CO2_error_quant(2:m, i);
                        best_mean_annual_precip = exp((bestCIAK+172.75)/35.7914); 
                        %This is the best fit line if MAP in regression where 
                        %MAP is considered as the independent variable. 
                        %If CIA-K is considered the independent variable, then 
                        %the best fit line is: 
                        %221.1*exp(0.02*paleosol_CO2_error_quant(2:m, i)) Sheldon et al 2002
                    end
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

                unexplained_sum_of_squares = sum((CIAK_Sheldon_calibration-yhat).^2);  %pg 469, also called error sum of squares

                s2YX = unexplained_sum_of_squares/(length(CIAK_Sheldon_calibration)-2);   %pg 473, same as MSE above from Davis. 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
    
        
        %paleoaltitude and paleolatitude needed for Kohn's 2010 equation to calculate
        %the big delta between atmosphere and plant leaves from MAP
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
    
    
        best_big_delta_leaf_atm = 2.01 - 1.98 *10^-4 * altitude + 5.88 * log10(best_mean_annual_precip + 300) + 0.00192 * abs(latitude); 
        %calcuates big delta (leaf-atm) using Kohn 2010
        bestd13Cp = (bestd13Ca - best_big_delta_leaf_atm)./(best_big_delta_leaf_atm./1000+1);
        %calculates d13C leaves from d13Ca (calculated above) and Kohn 2010 big delta
    
        bestd13Cr = (bestd13Cp - best_big_delta_leaf_srCO2)./(best_big_delta_leaf_srCO2./1000+1);  
        %Uses Bowling et al to determine d13Cr from d13C value of leaves
    
    
        %d13Cr error is calculated in MC subroutine for d13Crmethod = 4
    
    else  %if no MAP value available then can still estimate d13Cplant from d13Ca, but uncertainty is larger
        
        best_big_delta_leaf_atm = 19.1; %from Kohn, mean value for MAP < 1000 mm
        bestd13Cp = (bestd13Ca - best_big_delta_leaf_atm)./(best_big_delta_leaf_atm./1000+1);
        %calculates d13C leaves from d13Ca (calculated above) and Kohn 2010 big delta
        
        bestd13Cr = (bestd13Cp - best_big_delta_leaf_srCO2)./(best_big_delta_leaf_srCO2./1000+1);  
        %Uses Bowling et al to determine d13Cr from d13C value of leaves
    end

end