%PBUQ

%This program calculates the uncertainty associated with individual determinations of
%atmospheric CO2 made using Cerling's (1999) paleosol carbonate CO2 paleobarometer.
%The program considers error associated with measurements and also error
%associated with regression equations (e.g. transfer functions). A Monte Carlo simulation 
%is used to propagate uncertainty through the paleobarometer equation. 
%The program plots medians and percentile based estimates of confidence 
%intervals for atmospheric CO2.



%Data must be loaded as a text file named 'paleosol_CO2_error_quant.txt'.
%The file should have a separate column for each variable measured (it is
%not necessary to include values for every variable below in the text file).
%The first row of each column should contain a number (1-24) that indicates
%which variable is listed in the column, as follows
%1  for Age (Ma)
%2  for d13C value of paleosol carbonate
%3  for error on d13C value of paleosol carbonate
%4  for D47 temperature of carbonate formation
%5  for alkali ratio (Na2O + K2O)/Al2O3
%6  for other estimate of temperature
%7  for error associated with D47 or other estimate of temperature (the
        %error associated with alkali ratio dervied temperatures is assumed
        %constant (Retallack 2009)


%8  for depth to Bk (cm)
%9  for burial depth of paleosol (km)
%10 for CIA-K
%11 for other estimate of MAP (mm)
%12 for 1 sigma error on other estimate of MAP (mm)
%13 for other estimate of S(z)
%14 for 1 sigma error associated with other estimate of S(z)

%15 for d13C value of bulk paleosol organic matter
%16 for d13C value of organic matter occluded in paleosol carbonate or 1otherwise well preserved in paleosol
%17 for d13C value of contemporaneous coal or other minorly decomposed OM
        %from humid environment (for use with Arens et al (2000)
%18 for other estimate of d13C value of respired CO2
%19 for 1 sigma error associated with d13C value of OM or respired CO2
%20 for error associated with d13C value of coal or other minorly decomposed, humid environment OM

%21 for d13C value of contemporaneous marine carbonate
%22 for other estimate of d13Ca
%23 for 1 sigma error associated with d13C values of planktonic foram carbonate or other estimate of d13Ca
%24 for latitude (in degrees)
%25 for altitude (in m)

%26 for soil order (1 = Mollisol, 2 = Alfisol, 3 = Aridisol, 4 = Vertisol, 5 = Andisol, 6 = Inceptisol)
%27 for horizon SOM sampled from (1 = A, 2 = B or C)

clear all


%load data file with data from paleosols of interest
Table = readtable ('paleosol_CO2_error_quant.xlsx');
paleosol_CO2_error_quant = Table{:,:};

[m,n] = size(paleosol_CO2_error_quant);

%This loop pulls ages out of data file and creates an array
for i=1:n
    if paleosol_CO2_error_quant(1,i) == 1 %finds column with ages
        age = paleosol_CO2_error_quant(2:m,i);
    end 
end



%calculate the best values and define errors for the variables that go into the paleosol
%barometer. This is done with a separate subroutine for each variable in
%paleobarometer
delta13C_soilCO2
soil_derived_component_of_soil_CO2
delta13C_soil_respired_CO2
delta13C_atmospheric_CO2


%evaluate the expression for the paleosol barometer
bestCa = bestSz.*(bestd13Cs-1.0044*bestd13Cr-4.4)./(bestd13Ca-bestd13Cs);


%Monte carlo uncertainty simulation
Monte_Carlo_error_propagation



%Plot results
%First sort atm CO2 concentrations
sortedCa = sort(Ca);

%and then find percentile values
%middle 68 is comparable to 1 standard deviation for normal distributions
upper84percentCL_Ca=sortedCa(round(0.84*length(Ca)),:)';  %84th percntile value
lower16percentCL_Ca=sortedCa(round(0.16*length(Ca)),:)';  %16th percentile values

figure(5)
subplot(2,1,1)
errorbar(age, atm_CO2_estimate(:,1), atm_CO2_estimate(:,1)-lower16percentCL_Ca, upper84percentCL_Ca-atm_CO2_estimate(:,1), 'o-k');
hold on
xlabel('Age (Ma)')
ylabel('atmospheric CO2 concentration (ppmV)')

sortedRATIO = sort(SminusR_divby_AminusS);
upper84percentCL_RATIO=sortedRATIO(round(0.84*length(SminusR_divby_AminusS)),:)';
lower16percentCL_RATIO=sortedRATIO(round(0.16*length(SminusR_divby_AminusS)),:)';

subplot(2,1,2)
errorbar(age, RATIO, RATIO-lower16percentCL_RATIO, upper84percentCL_RATIO-RATIO, 'o-b');
hold on
xlabel('Age (Ma)')
ylabel('RATIO')
