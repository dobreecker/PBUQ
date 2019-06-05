%This is a subroutine within PBUQ

%The soil-derived component of total soil CO2, S(z), is associated with
%larger relative errors than any other the other variables used to
%determine atm CO2 cocnentration from paleosol carbonate. Cotton and
%Sheldon 2012 have developed a relationship betwwen mean annual
%precipitation (MAP) and S(z). This program allows users to use the Cotton and
%Sheldon 2012 calibration (and consider its associated error) if MAP can be
%determined somehow. The program incorporates several option for
%determining MAP. Montanez (2013) has S(z) defined ranges for various paleosol
%types- those ranges have been modified and are another option for
%estimated S(z).


%ask user to select the method preferred for S(z)
Szmethod = input ('How do you want to determine S(z)? \n 1) from MAP (Cotton and Sheldon 2012) \n 2) from depth to Bk (Retallack 2009) \n 3) S(z) based on soil order (Montanez 2013) \n 4) other estimate of S(z) \n')

if Szmethod == 1 %S(z) from MAP according to Cotton and Sheldon 2012. 
                 %If user wants to calculate S(z) from MAP, then need to 
                 %know how user wants to determine MAP
    
    MAPmethod = input ('How do you want to determine MAP for calculating S(z)? \n 1) from depth to B(k) (Retallack 2005) \n 2) from CIA-K (Sheldon et al 2002) \n 3) other estimate of MAP \n')
    
    if MAPmethod == 1 %MAP from depth to Bk. See subroutine 'd13C_soilCO2' 
                      %for a description of what loops like this do
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
        %considered here, but when propogated are likely smaller than other
        %errors
        best_decompacted_depth_to_Bk = best_depth_to_Bk./(-0.62./(0.38./exp(0.17*best_burial_depth)-1));
        %equation for decompaction for Aridisols (Sheldon and Tabor 2009)
        
        log_best_decompacted_depth_to_Bk = log(best_decompacted_depth_to_Bk);
        %log transformation used before fitting regression line
        
        best_mean_annual_precip = (log_best_decompacted_depth_to_Bk-2.97)/.0021; 
        %This is the best fit regression line if MAP is considered the 
        %independent variable.  If depth to Bk considered the independent 
        %variable then the best fit polynomial is 
        %137.24 + 6.45*best_decompacted_depth_to_Bk - 
        %0.013*best_decompacted_depth_to_Bk.^2 (Retallack 2005)
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %***********X FROM Y**********
        %Quantifying the error associated with estimating X from Y (in this case
        %estimating MAP from depth to Bk) is discussed by
        %Sokal and Rohlf 1981 'Biometry'. PBUQ uses an expression for 
        %confidence limits to estimate pdfs for X. This calculation is 
        %performed in the Monte Carlo subroutine,but the required preliminary 
        %calculations are carried out here because these don't require new 
        %calculations for each iteration

        load Retallack_2005_data.txt     %load calibration data

        MAP_R_2005 = Retallack_2005_data(:,1);
        depth_to_Bk_R_2005 = Retallack_2005_data(:,2);
        log_d_to_Bk = log(depth_to_Bk_R_2005); %natural log (log transformation of depth to Bk for regression)

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
                best_mean_annual_precip = exp((bestCIAK+172.75)/35.7914); 
                %This is the best fit line in regression where MAP is 
                %considered as the independent variable. If CIA-K is considered 
                %the independent variable, then the best fit line is: 
                %221.1*exp(0.02*paleosol_CO2_error_quant(2:m, i)) Sheldon et al 2002
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %***********X FROM Y (see lines 57-63 for explanation)**********

        load Sheldonetal_2002_CIAKdata.txt      %load calibration data

        ln_MAP_Sheldon_calibration = Sheldonetal_2002_CIAKdata(:,1);
        CIAK_Sheldon_calibration = Sheldonetal_2002_CIAKdata(:,2);


        %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
        sum_of_squares_of_Y = sum((CIAK_Sheldon_calibration-mean(CIAK_Sheldon_calibration)).^2); %pg 55
        sum_of_squares_of_X = sum((ln_MAP_Sheldon_calibration-mean(ln_MAP_Sheldon_calibration)).^2);     %pg 55
        sum_xy = sum((ln_MAP_Sheldon_calibration-mean(ln_MAP_Sheldon_calibration)).*...
            (CIAK_Sheldon_calibration-mean(CIAK_Sheldon_calibration)));  %Eq 14.4 page 467

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
    
    
    
    %Now that the best value for MAP has been defined, use Cotton and
    %Sheldon 2012 to calculate a best value for S(z)
    bestSz = (5.67*best_mean_annual_precip-269.9);   
    
    %And now determine error associated with S(z)
    Sz_error_calculator 
    %This subroutine calculates error associated with Cotton and Sheldon 2012 regression line
    
    SE_Sztransferfunctionerror = sqrt(MSE*(1 + 1/length_Cotton_SheldonSz + (best_mean_annual_precip-meanCandSMAP).^2/sum_xi_minus_xbar_squared));
        %Error for regression type 'Y from X.'
        %This uses values calculated in subroutine Sz_error_calculator with
        %values calculated above for an individual paleosol to calculate
        %standard error at the appropriate segment (for the paleosol of interest) of the transfer function 
        %This is the standard error of a new observation (not standard
        %error of the mean).

elseif Szmethod == 2  %S(z) directly from depth to Bk (Retallack 2009)
    
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
    %considered here, but when propogated are likely smaller than other
    %errors
    best_decompacted_depth_to_Bk = best_depth_to_Bk./(-0.62./(0.38./exp(0.17*best_burial_depth)-1));
    
    %bestSz = 67.0*best_decompacted_depth_to_Bk+567;
    
    %revised to reflect minimum warm season S(z)
    bestSz = 35.3*best_decompacted_depth_to_Bk+588;
   
    %And now determine error associated with S(z)
    Sz_from_depth_to_Bk_error_calculator
    %This subroutine calculates error associated with Retallack 2009 regression line
    
    SE_Sz_from_Bk_transferfunctionerror = sqrt(MSE*(1 + 1/length_Sz_Retallack_2009 + (best_burial_depth-mean_Retallack_2009_depth_to_Bk).^2/sum_xi_minus_xbar_squared));
        %Error for regression type 'Y from X.'
        %This uses values calculated in subroutine Sz_from_depth_to_Bk_error_calculator with
        %values calculated above for an individual paleosol to calculate
        %standard error at the appropriate segment (for the paleosol of interest) of the transfer function 
        %This is the standard error of a new observation (not standard
        %error of the mean).
        
        
%original PBUQ code below- S(z) as independent variable
%     bestSz = (best_decompacted_depth_to_Bk - 1.632)/0.012;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %***********X FROM Y (see lines 57-63 for explanation)**********
%     
%     Sz_Retallack_2009 = [7921, 5616, 4342, 4474, 4274, 4614, 2374, 2853, 2274, 2417, 3938, 2339, 1404, 664, 1153];
%     depth_to_Bk_Retallack_2009 = [110, 78, 66, 40, 40, 40, 30, 36, 40, 28, 40, 20, 13, 13, 35];
% 
%     %Using notation in Sokal and Rohlf- This is for one value of Y for each value of X
%     sum_of_squares_of_Y = sum((depth_to_Bk_Retallack_2009-mean(depth_to_Bk_Retallack_2009)).^2); %pg 55
%     sum_of_squares_of_X = sum((Sz_Retallack_2009-mean(Sz_Retallack_2009)).^2);     %pg 55
%     sum_xy = sum((Sz_Retallack_2009-mean(Sz_Retallack_2009)).*(depth_to_Bk_Retallack_2009-mean(depth_to_Bk_Retallack_2009)));  %Eq 14.4 page 467
% 
%     b_yx = sum_xy/sum_of_squares_of_X; %slope of regression line
% 
%     Y_intercept = mean(depth_to_Bk_Retallack_2009)-b_yx*mean(Sz_Retallack_2009);
% 
%     yhat = b_yx * Sz_Retallack_2009 + Y_intercept;
% 
%     unexplained_sum_of_squares = sum((depth_to_Bk_Retallack_2009-yhat).^2);  %pg 469, also called error sum of squares
% 
%     s2YX = unexplained_sum_of_squares/(length(depth_to_Bk_Retallack_2009)-2);   %pg 473, same as MSE above from Davis. 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of original PBUQ removed 
    

elseif Szmethod == 3 %S(z) from soil orders (modified from Montanez 2013)
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 26   %The column which identifies the type of paleosol (by number)
            soilorder = paleosol_CO2_error_quant(2:m, i);
        end
    end
    
    %Modified S(z) by soil are as follows
    
    Mollisol_Sz = [365
    417
    517
    299
    247
    2446
    479
    611
    593
    398
    650
    493
    405
    846
    881
    631
    1132
    904
    627
    754
    735
    550
    517
    1059
    462
    424
    643
    609
    583
    450
    1101
    944
    615
    527
    644
    491
    449
    820
    596
    487
    501
    488
    465
    522
    599
    518
    518
    404
    399
    367
    415
    465
    447
    422
    539
    668
    550
    450
    475
    446
    465
    668
    370
    719
    763
    11577
    610
    349
    636
    1035
    1204
    971
    586
    520
    439
    482
    461
    737
    639
    895
    1556
    1019
    1249
    1001
    623
    735
    953
    1825
    1687
    996
    1022
    803
    1215];
    
    Alfisol_Sz = [22
    983
    544
    792
    681
    532
    641
    792
    602
    548
    538
    558
    499
    629
    797
    1064
    873];


    Aridisol_Sz = [357;
    315
    900
    1309
    458
    1183
    885
    344
    202
    964
    741
    555
    440
    355
    381
    237
    437
    357
    298
    283
    248
    277
    257
    694
    586
    618
    1549
    1305
    724
    2592
    1192
    606
    159
    63
    79];

    Vertisol_Sz = [1536
    547
    1270
    1465
    753
    1270
    1287
    711
    380
    503
    209
    4571
    2355
    546
    3474
    488
    229
    901
    540
    1081
    687
    634
    268
    191
    341
    3665
    7547
    1653
    712
    244
    224
    232
    2921
    280
    432
    339
    159
    7294
    2631
    2670
    4503];


    Anidsol_Sz = [2867
    2359
    1538
    1451
    1328
    806
    666
    540
    1211
    1291];
    
    Inceptisol_Sz = [327
    320
    293
    356
    153
    269
    149
    216
    515
    372
    496
    268
    270
    224
    282
    329];
    
    %hitograms of S(z) values are skewed right, so medians are used as best
    %S(z) values to avoid biasing the best S(z) with a few very high S(z)
    for j = 1:m-1
        if soilorder(j) == 1
            bestSz(j) = median(Mollisol_Sz);
        elseif soilorder(j) == 2
            bestSz(j) = median(Alfisol_Sz);
        elseif soilorder(j) == 3
            bestSz(j) = median(Aridisol_Sz);
        elseif soilorder(j) == 4
            bestSz(j) = median(Vertisol_Sz);
        elseif soilorder(j) == 5
            bestSz(j) = median(Andisol_Sz);
        elseif soilorder(j) == 6
            bestSz(j) = median(Inceptisol_Sz);
        end
        
    end
    bestSz = bestSz';
    
elseif Szmethod ==4
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 13
            bestSz = paleosol_CO2_error_quant(2:m, i);
        end
    end
    
    for i=1:n 
        if paleosol_CO2_error_quant(1,i) == 14
            SD_Sz = paleosol_CO2_error_quant(2:m, i);
        end
    end
    
end

