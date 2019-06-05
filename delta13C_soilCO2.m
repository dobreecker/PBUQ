%This is a subroutine within PBUQ

%The d13C value of ancient soil CO2 is calculated from the d13C value of 
%paleosol carbonate and the temperature of carbonate formation. Carbon isotope
%equilibrium between soil CO2 and soil carbonate is assumed. This
%program asks user which method is preferred to determine the temperature
%of soil carbonate formation. D47 values are directly related to the temperature of
%carbonate formation and should be used when available. Estimates of mean
%annual air temperature can also be used to determine temperature of soil
%carbonate formation, using modern calibrations reported by Quade and others
%and Passey and others (2010).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First determine temperature and its associated error


%ask user to select the method preferred for temperature
tempcalcmethod = input ('How do you want to determine temperature of carbonate formation? \n 1) D47 \n 2) alkali ratio \n 3) other mean annual air temperature estimate \n')


if tempcalcmethod == 1
    %search through the first row of text file until find '4', which is the 
    %column heading for D47 temperatures
    
    for i=1:n %n is number of columns in the input text file
        if paleosol_CO2_error_quant(1,i) == 4
            besttemperature = paleosol_CO2_error_quant(2:m, i);   
            %row 2 through the last row (row m) of the column with heading 
            %'4' contains best estimates of temperature (row 1 has the 
            %column label, in this case the number '4')
        end
    end
    MAT = (besttemperature-17.974)./0.506;   
    %from compilation of Quade et al in press and Passey et al 2010 
    %comparisons of soil carbonate D47 T and MAAT             
    
    %MAT (mean annual air temperature) is calculated here because then the 
    %same Monte Carlo subroutine (which runs next and uses MAT values 
    %calculated here) can be run for all options (tempcalcmethod 1-3). 
    %The Monte Carlo subroutine calculates carbonate formation T from
    %MAT. For this option (tempcalcmethod = 1) in which MAT is determined 
    %from D47 measurements, no error should associated with
    %making this calculation again in reverse (i.e. calculation of
    %carbonate formation T from MAT), so we write the following:
    SE_MAATtocarbformTtransferfunctionerror = 0;  
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 7 
            %column 7 contains standard deviations associated with D47 
            %temperature determinations
            SD_MAT = paleosol_CO2_error_quant(2:m, i);
        end
    end   
    %The SD_MAT values above are really the standard error of carbonate 
    %formation temperatures (i.e. D47 temperatures), not MAT, but since 
    %MATs are being put in to the Monte Carlo subroutine, and no error will 
    %be added to this (SE_MAATtocarbformTtransferfunctionerror = 0), it is 
    %convenient to input this as SD_MAT, which is what the Monte Carlo 
    %subroutine expects for the next two tempcalcmethod options
    
    
elseif tempcalcmethod == 2   %for alkali ratio thermometer
    for i=1:n
        if paleosol_CO2_error_quant(1,i) == 5
            MAT = -18.516*paleosol_CO2_error_quant(2:m, i)+17.298;       %Sheldon et al 2002
        end
    end
    besttemperature = 0.506 * MAT + 17.974;      
    %from compilation of Quade et al in review and Passey et al 2010 
    %comparisons of soil carb D47 T and MAAT
    
    SD_MAT = 4.4;                                
    %I have not quantified this using an error calculator subroutine, like 
    %for other transfer functions, instead I use the value reported by 
    %Sheldon et al 2002
    
    carbformT_error_calculator 
    %This subroutine calculates the error associated with the regression 
    %equation relating mean annual air temeprature and carbonate formation 
    %temperature
    
    SE_MAATtocarbformTtransferfunctionerror = sqrt(MSE*(1 + 1/lengthD47T +...
        (MAT-mean_mean_annual_air_temperature).^2/sum_xi_minus_xbar_squared));
    %Error for regression type 'Y from X.' 
    %MAT is mean annual temperature determined for the paleosol(s) of interest,
    %mean_mean_annual_air_temperature is the mean of MAATs investigated by
    %Quade et al 2013 and Passey et al 2010 in their calibration studies 
    %of the relationship between MAAT and soil carbonate formation temperature.
    %This equation is for the standard error of new onservations 
    %(NOT the standard error of the mean) from Davis 2002 
    %'Statistics and Data Analysis in Geology' pg. 204.

   
elseif tempcalcmethod == 3  %for other user-specified MAT estimate. This is if user wants to use a different proxy
    for i=1:n
        if paleosol_CO2_error_quant(1,i) == 6
            MAT = paleosol_CO2_error_quant(2:m, i);
        end
    end
    besttemperature = 0.506 * MAT + 17.974;      
    %from compilation of Quade et al in review and Passey et al 2010 
    %comparisons of soil carb D47 T and MAAT
    
    carbformT_error_calculator
    %This subroutine calculates the error associated with the regression 
    %equation relating mean annual air temeprature and carbonate formation 
    %temperature
    
    SE_MAATtocarbformTtransferfunctionerror = sqrt(MSE*(1 + 1/lengthD47T + (MAT-mean_mean_annual_air_temperature).^2/sum_xi_minus_xbar_squared));
    %Error for regression type 'Y from X.'
    %MAT is mean annual temperature determined for the paleosol(s) of interest,
    %mean_mean_annual_air_temperature is the mean of MAATs investigated by
    %Quade et al 2013 and Passey et al 2010 in their calibration studies 
    %of the relationship between MAAT and soil carbonate formation temperature.
    %This equation is for the standard error of new onservations 
    %(NOT the standard error of the mean) from Davis 2002 
    %'Statistics and Data Analysis in Geology' pg. 204.
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 7 
            %column 7 contains standard deviations associated with temperature 
            %(for this option the errors listed in text file should be 
            %associated with estimates of MAAT estimated using some other 
            %proxy not considered by PBUQ)
            SD_MAT = paleosol_CO2_error_quant(2:m, i);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now define the d13C value of paleosol carbonate and its associated error


for i=1:n
    if paleosol_CO2_error_quant(1,i) == 2 %finds row with d13C values of paleosol carbonate   
        bestd13Cpc = paleosol_CO2_error_quant(2:m, i);
    end
end
  
for i=1:n
    if paleosol_CO2_error_quant(1,i) == 3 %finds row with standard deviations of measured d13C values of paleosol carbonate
        SD_d13Cc = paleosol_CO2_error_quant(2:m, i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Finally, calculate d13Cs using the values determined above
 bestd13Cs = (bestd13Cpc + 1000)./((11.98-0.12*besttemperature)./1000+1)-1000;   
%Romanek et at 1992- error associated with constants in this equation not 
%considered here, but likely to be much smaller than other errors
