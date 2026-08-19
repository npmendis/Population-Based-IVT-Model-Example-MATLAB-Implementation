%% 
% Load lookup table for ionic equilibria
load('gridded_lookup_table_5_50_50_50.mat', 'X1', 'X2', 'X3', 'X4', 'output_data_gridded');
F_interp = griddedInterpolant(X1, X2, X3, X4, output_data_gridded, 'linear');

% Define constant initial conditions
C_T7RNAP_init = [0.4e-3, 0.4e-3, 0.4e-3]; % mM
C_DNAP_init = [0.5e-3, 0.5e-3, 0.5e-3]; % mM
C_ATP_total_init = [0, 0, 0]; % mM
C_UTP_total_init = [4, 6, 6]; % mM
C_CTP_total_init = [4, 6, 6]; % mM
C_GTP_total_init = [4, 6, 6]; % mM
C_Mg_total_init = [20, 20, 30]; % mM
C_H_total_init = [15.4812, 15.5104, 15.4825]; % mM (calculated)

% experimental data from literature
t_rna12 = [3.72;
           7.06;
           15.3;
           32.1;
           64.1;
           128;
           180;
           241
           ] / 60; %in h

rna12 =   [0.0583;
           0.154;
           0.24;
           0.341;
           0.337;
           0.317;
           0.312;
           0.314
          ]; %in mM

t_abt12 = [2.97;
           7.2;
           15.9;
           31.7;
           64.1;
           128;
           180;
           240
          ] / 60; %in h

abt12 =   [0.154;
           0.363;
           0.578;
           0.885;
           0.968;
           1;
           0.999;
           1.02
          ]; %in mM

t_rna12_6 = [8.24;
             15.9;
             32.4;
             64.7;
             122;
             180;
             241
           ] / 60; %in h

rna12_6 =   [0.0938;
             0.204;
             0.322;
             0.307;
             0.317;
             0.246;
             0.363
          ]; %in mM

t_abt12_6 = [16.1;
             32.9;
             64.9;
             121;
             180;
             241
          ] / 60; %in h

abt12_6 =   [0.784;
             1.21;
             1.25;
             1.34;
             1.54;
             1.5
          ]; %in mM

t_rna12_30 = [7.86;
             16.3;
             32.2;
             64.7;
             129;
             180;
             241
           ] / 60; %in h

rna12_30 =   [0.126;
              0.345;
              0.397;
              0.538;
              0.502;
              0.517;
              0.457
          ]; %in mM

t_abt12_30 = [8.18;
              16.7;
              32.4;
              64.9;
              129;
              181;
              241            
          ] / 60; %in h

abt12_30 =   [0.301;
              0.724;
              0.96;
              1.26;
              1.38;
              1.34;
              1.44
          ]; %in mM

% Number of experiments
num_experiments = 1;

% Store experimental data in a structured format
exp_data = struct();
for i = 1:num_experiments

    exp_data(i).t_RNA12 = t_rna12(:,i);
    exp_data(i).RNA12 = rna12(:,i);
    exp_data(i).t_abt12 = t_abt12(:,i);
    exp_data(i).abt12 = abt12(:,i);
    exp_data(i).t_RNA12_6 = t_rna12_6(:,i);
    exp_data(i).RNA12_6 = rna12_6(:,i);
    exp_data(i).t_abt12_6 = t_abt12_6(:,i);
    exp_data(i).abt12_6 = abt12_6(:,i);
    exp_data(i).t_RNA12_30 = t_rna12_30(:,i);
    exp_data(i).RNA12_30 = rna12_30(:,i);
    exp_data(i).t_abt12_30 = t_abt12_30(:,i);
    exp_data(i).abt12_30 = abt12_30(:,i);

end

%optimal parameter values in the following order:
% [k_i, ke, k_a, K_1, K_2, k_on, k_off]
best_params = [4.756 8.945 7.855 0.3601 -0.1672 7.477 2.046];

% ---------- PLOTTING RESULTS ----------

%12mer

figure;
for i = 1:num_experiments
    % Solve ODE with estimated parameters
    tspan = linspace(0,16200/3600,100); % All unique time points
    [t_model_12, C_model_12] = solve_transcription_population_model_12mer(best_params, ...
        C_T7RNAP_init(1), C_DNAP_init(1), C_ATP_total_init(1), ...
        C_UTP_total_init(1), C_CTP_total_init(1), C_GTP_total_init(1), ...
        C_Mg_total_init(1), C_H_total_init(1), tspan, F_interp);

    [t_model_12_6, C_model_12_6] = solve_transcription_population_model_12mer(best_params, ...
        C_T7RNAP_init(2), C_DNAP_init(2), C_ATP_total_init(2), ...
        C_UTP_total_init(2), C_CTP_total_init(2), C_GTP_total_init(2), ...
        C_Mg_total_init(2), C_H_total_init(2), tspan, F_interp);

    [t_model_12_30, C_model_12_30] = solve_transcription_population_model_12mer(best_params, ...
        C_T7RNAP_init(3), C_DNAP_init(3), C_ATP_total_init(3), ...
        C_UTP_total_init(3), C_CTP_total_init(3), C_GTP_total_init(3), ...
        C_Mg_total_init(3), C_H_total_init(3), tspan, F_interp);

    % Extract model predictions for RNA and PPi
    RNA_model_12 = C_model_12(:, 16) + C_model_12(:, 15);
    abt_model_12 = C_model_12(:, 17);

    for j=5:1:12
        abt_model_12 = abt_model_12 + C_model_12(:, j);
    end

    RNA_model_12_6 = C_model_12_6(:, 16) + C_model_12_6(:, 15);
    abt_model_12_6 = C_model_12_6(:, 17);

    for j=5:1:12
        abt_model_12_6 = abt_model_12_6 + C_model_12_6(:, j);
    end

    RNA_model_12_30 = C_model_12_30(:, 16) + C_model_12_30(:, 15);
    abt_model_12_30 = C_model_12_30(:, 17);

    for j=5:1:12
        abt_model_12_30 = abt_model_12_30 + C_model_12_30(:, j);
    end

    % Interpolate model predictions to experimental time points
    RNA_model_12_interp = interp1(t_model_12, RNA_model_12, exp_data(i).t_RNA12, 'linear');
    abt_model_12_interp = interp1(t_model_12, abt_model_12, exp_data(i).t_abt12, 'linear');

    RNA_model_12_6_interp = interp1(t_model_12_6, RNA_model_12_6, exp_data(i).t_RNA12_6, 'linear');
    abt_model_12_6_interp = interp1(t_model_12_6, abt_model_12_6, exp_data(i).t_abt12_6, 'linear');

    RNA_model_12_30_interp = interp1(t_model_12_30, RNA_model_12_30, exp_data(i).t_RNA12_30, 'linear');
    abt_model_12_30_interp = interp1(t_model_12_30, abt_model_12_30, exp_data(i).t_abt12_30, 'linear');

    RNA_exp_all = [exp_data(i).RNA12; exp_data(i).RNA12_6; exp_data(i).RNA12_30; exp_data(i).abt12; exp_data(i).abt12_6; exp_data(i).abt12_30];
    RNA_model_interp_all = [RNA_model_12_interp; RNA_model_12_6_interp; RNA_model_12_30_interp; abt_model_12_interp; abt_model_12_6_interp; abt_model_12_30_interp];

    % nAAD (%)
    nAAD_J1 = 100*mean(abs(RNA_model_interp_all - RNA_exp_all) / max(RNA_exp_all));

    fprintf('Overall J1 nAAD = %.2f %%\n', nAAD_J1);

    % Combined RNA and PPi comparison
    subplot(1, num_experiments, i)
    hold on;

    % Plot experimental RNA: blue filled diamond
    plot(exp_data(i).t_RNA12, exp_data(i).RNA12, ...
        's', 'MarkerSize', 8, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '[0 0.25 0.90]', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Runoff RNA (Mg: 20 mM, NTP: 4 mM each)');

    plot(exp_data(i).t_RNA12_6, exp_data(i).RNA12_6, ...
        's', 'MarkerSize', 8, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '[1.00 0.00 0.00]', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Runoff RNA (Mg: 20 mM, NTP: 6 mM each)');

    plot(exp_data(i).t_RNA12_30, exp_data(i).RNA12_30, ...
        's', 'MarkerSize', 8, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', '[0.00 0.50 0.00]', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Runoff RNA (Mg: 30 mM, NTP: 6 mM each)');
    
    % Plot model RNA: solid orange line
    plot(t_model_12, RNA_model_12, ...
        '-', 'LineWidth', 2.7, 'Color', '[0 0.25 0.90]', 'HandleVisibility', 'off');

    plot(t_model_12_6, RNA_model_12_6, ...
        '-', 'LineWidth', 2.7, 'Color', '[1.00 0.00 0.00]', 'HandleVisibility', 'off');    

    plot(t_model_12_30, RNA_model_12_30, ...
        '-', 'LineWidth', 2.7, 'Color', '[0.00 0.50 0.00]', 'HandleVisibility', 'off');    

    
    % Plot experimental aborts: hollow blue square
    plot(exp_data(i).t_abt12, exp_data(i).abt12, ...
        'd', 'MarkerSize', 8, 'MarkerEdgeColor', '[0 0.25 0.90]', 'MarkerFaceColor', 'none', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Aborts (Mg: 20 mM, NTP: 4 mM each)');

    plot(exp_data(i).t_abt12_6, exp_data(i).abt12_6, ...
        'd', 'MarkerSize', 8, 'MarkerEdgeColor', '[1.00 0.00 0.00]', 'MarkerFaceColor', 'none', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Aborts (Mg: 20 mM, NTP: 6 mM each)');    

     plot(exp_data(i).t_abt12_30, exp_data(i).abt12_30, ...
        'd', 'MarkerSize', 8, 'MarkerEdgeColor', '[0.00 0.50 0.00]', 'MarkerFaceColor', 'none', 'LineStyle', 'none', 'LineWidth', 1.4, ...
        'DisplayName', 'Aborts (Mg: 30 mM, NTP: 6 mM each)');   

    % Plot model aborts: dashed orange line
    plot(t_model_12, abt_model_12, ...
        '--', 'LineWidth', 2.7, 'Color', '[0 0.25 0.90]', 'HandleVisibility', 'off'); 

    plot(t_model_12_6, abt_model_12_6, ...
        '--', 'LineWidth', 2.7, 'Color', '[1.00 0.00 0.00]', 'HandleVisibility', 'off'); 

    plot(t_model_12_30, abt_model_12_30, ...
        '--', 'LineWidth', 2.7, 'Color', '[0.00 0.50 0.00]', 'HandleVisibility', 'off');     

    xlabel('Time (h)');
    ylabel('Concentration (mM)');
    %title(sprintf('Exp %d: RNA & PPi', i));
    %ylim([0 50]);

    legend;
    lgd = legend('Location', 'best');
    lgd.FontSize = 7;
    lgd.Box = 'off';

    hold off;   
end

function [t, C] = solve_transcription_population_model_12mer(params, ...
    C_T7RNAP_init, C_DNAP_init, C_ATP_total_init, C_UTP_total_init, C_CTP_total_init, C_GTP_total_init, C_Mg_total_init, C_H_total_init,...
    tspan_init, F_interp)

    %define RNA sequence
    RNA = 'GGCGCUUGCGUC'; % G5 A0 C4 U3

    k_on = 10^params(6); % h-1mM-1
    k_off = 10^params(7); % h-1
    k_i = 10^params(1); % h-1
    k_a = 10^params(3); %h-1
    k_e1 = 10^params(2); % h-1
    k_e2 = 10^params(2); % h-1
    K_1 = 10^params(4); % mM
    K_2 = 10^params(5); % mM

    % Initial conditions
    C_T7RNAP_0 = C_T7RNAP_init; %mM
    C_T7RNAP_DNAP_0 = 0; %mM
    C_DNAP_0 = C_DNAP_init; %mM
    C_T7RNAP_M_0 = zeros(1,12); %mM
    C_RNA_0 = 0; %mM
    C_imp_0 = 0; %mM
    C_ATP_total_0 = C_ATP_total_init; %mM
    C_UTP_total_0 = C_UTP_total_init; %mM
    C_CTP_total_0 = C_CTP_total_init; %mM
    C_GTP_total_0 = C_GTP_total_init; %mM
    C_Mg_total_0 = C_Mg_total_init; %mM
    C_PPi_total_0 = 0; %mM
    C_H_total_0 = C_H_total_init; %mM
    
    % Time span
    tspan = tspan_init;

    % Solve ODEs
    [t, C] = ode15s(@(t, C) transcription_population_model_12mer(t, C, k_on, k_off, k_i, k_a, k_e1, k_e2, K_1, K_2, F_interp),...
        tspan, [C_T7RNAP_0, C_T7RNAP_DNAP_0, C_DNAP_0, C_T7RNAP_M_0, C_RNA_0, C_imp_0, C_ATP_total_0, C_UTP_total_0, C_CTP_total_0, C_GTP_total_0, C_Mg_total_0, C_PPi_total_0, C_H_total_0]);

end

function dCdt = transcription_population_model_12mer(t, C, k_on, k_off, k_i, k_a, k_e1, k_e2, K_1, K_2, F_interp)
    
    %define RNA sequence
    RNA = 'GGCGCUUGCGUC'; % G5 A0 C4 U3 

    %unpack variables
    C_T7RNAP = C(1);
    C_T7RNAP_DNAP = C(2); 
    C_DNAP = C(3);
       
    C_T7RNAP_M = zeros(12,1);
    C_T7RNAP_M = C(4:15);
    
    C_RNA = C(16);
    C_imp = C(17);

    C_ATP_total = C(18);
    C_UTP_total = C(19);
    C_CTP_total = C(20);
    C_GTP_total = C(21);

    C_Mg_total = C(22);

    C_PPi_total = C(23);
    C_H_total = C(24);    

    C_NTP_total = C_ATP_total + C_UTP_total + C_CTP_total + C_GTP_total;

    % Solve using lookup table
    solution = F_interp(C_Mg_total, C_NTP_total, C_PPi_total, C_H_total);    

    C_MgATP = solution(2) * C_ATP_total / C_NTP_total; % mM Example assumed constant (modify as needed)
    C_MgUTP = solution(2) * C_UTP_total / C_NTP_total; % mM Example assumed constant (modify as needed)
    C_MgCTP = solution(2) * C_CTP_total / C_NTP_total; % mM Example assumed constant (modify as needed)
    C_MgGTP = solution(2) * C_GTP_total / C_NTP_total ;% mM Example assumed constant (modify as needed)
    C_MgPPi = solution(5);
    C_Mg = solution(1);
    C_Mg2PPi = solution(6);

    %assign NTPs
    for i = 1:12
        if RNA(i) == 'A'
            NTP(i) = C_MgATP;
            ATP_on(i) = 1;
            UTP_on(i) = 0;
            CTP_on(i) = 0;
            GTP_on(i) = 0;
        elseif RNA(i) == 'U'
            NTP(i) = C_MgUTP;
            ATP_on(i) = 0;
            UTP_on(i) = 1;
            CTP_on(i) = 0;
            GTP_on(i) = 0;            
        elseif RNA(i) == 'C'
            NTP(i) = C_MgCTP;
            ATP_on(i) = 0;
            UTP_on(i) = 0;
            CTP_on(i) = 1;
            GTP_on(i) = 0;            
        else
            NTP(i) = C_MgGTP;
            ATP_on(i) = 0;
            UTP_on(i) = 0;
            CTP_on(i) = 0;
            GTP_on(i) = 1;            
        end
    end
    
    %define Michaelis-Menten terms
    f_MM(1) = NTP(1) * NTP(2) / (K_1^2 + K_1 * NTP(1) + K_1 * NTP(2) + NTP(1) * NTP(2)) * C_Mg / (K_2 + C_Mg);

    for i = 3:12
            f_MM(i) = NTP(i) / (K_1 + NTP(i)) * C_Mg / (K_2 + C_Mg);
    end

    % define ODEs
    dCdt = zeros(24,1);

    dCdt(1) = - k_on * C_T7RNAP * C_DNAP + k_off * C_T7RNAP_DNAP + k_e2 * C_T7RNAP_M(12); %E
    
    dCdt(2) = k_on * C_T7RNAP * C_DNAP - k_off * C_T7RNAP_DNAP - k_i * C_T7RNAP_DNAP * f_MM(1) + sum(k_a * C_T7RNAP_M(2:9)); %E-DNA
    
    dCdt(3) = - k_on * C_T7RNAP * C_DNAP + k_off * C_T7RNAP_DNAP + k_e2 * C_T7RNAP_M(12); %DNA

    dCdt(4) = 0; %E-DNA-T1

    dCdt(5) = k_i * C_T7RNAP_DNAP * f_MM(1) - k_e1 * C_T7RNAP_M(2) * f_MM(3) - k_a * C_T7RNAP_M(2); %%E-DNA-T2

    for i = 3:9
        dCdt(i+3) = k_e1 * C_T7RNAP_M(i-1) * f_MM(i) - k_e1 * C_T7RNAP_M(i) * f_MM(i+1) - k_a * C_T7RNAP_M(i); %%E-DNA-Ti
    end

    dCdt(13) = k_e1 * C_T7RNAP_M(9) * f_MM(10) - k_e2 * C_T7RNAP_M(10) * f_MM(11); %%E-DNA-T10

    dCdt(14) = k_e2 * C_T7RNAP_M(10) * f_MM(11) - k_e2 * C_T7RNAP_M(11) * f_MM(12); %%E-DNA-T11

    dCdt(15) = k_e2 * C_T7RNAP_M(11) * f_MM(12) - k_e2 * C_T7RNAP_M(12); %%E-DNA-T12

    dCdt(16) = k_e2 * C_T7RNAP_M(12); %RNA balance

    dCdt(17) = sum(k_a * C_T7RNAP_M(2:9)); %total abortive products

    dCdt(18) = - k_i * C_T7RNAP_DNAP * f_MM(1) * ATP_on(1) - k_i * C_T7RNAP_DNAP * f_MM(1) * ATP_on(2); %ATP
    dCdt(19) = - k_i * C_T7RNAP_DNAP * f_MM(1) * UTP_on(1) - k_i * C_T7RNAP_DNAP * f_MM(1) * UTP_on(2); %UTP
    dCdt(20) = - k_i * C_T7RNAP_DNAP * f_MM(1) * CTP_on(1) - k_i * C_T7RNAP_DNAP * f_MM(1) * CTP_on(2); %CTP
    dCdt(21) = - k_i * C_T7RNAP_DNAP * f_MM(1) * GTP_on(1) - k_i * C_T7RNAP_DNAP * f_MM(1) * GTP_on(2); %GTP
    dCdt(23) = k_i * C_T7RNAP_DNAP * f_MM(1); %PPi
    dCdt(24) = k_i * C_T7RNAP_DNAP * f_MM(1); %H

    for j = 3:10
        dCdt(18) = dCdt(18) - k_e1 * C_T7RNAP_M(j-1) * f_MM(j) * ATP_on(j); %ATP
        dCdt(19) = dCdt(19) - k_e1 * C_T7RNAP_M(j-1) * f_MM(j) * UTP_on(j); %UTP
        dCdt(20) = dCdt(20) - k_e1 * C_T7RNAP_M(j-1) * f_MM(j) * CTP_on(j); %CTP
        dCdt(21) = dCdt(21) - k_e1 * C_T7RNAP_M(j-1) * f_MM(j) * GTP_on(j); %GTP
        dCdt(23) = dCdt(23) + k_e1 * C_T7RNAP_M(j-1) * f_MM(j); %PPi
        dCdt(24) = dCdt(24) + k_e1 * C_T7RNAP_M(j-1) * f_MM(j); %H
    end

    for j = 11:12
        dCdt(18) = dCdt(18) - k_e2 * C_T7RNAP_M(j-1) * f_MM(j) * ATP_on(j); %ATP
        dCdt(19) = dCdt(19) - k_e2 * C_T7RNAP_M(j-1) * f_MM(j) * UTP_on(j); %UTP
        dCdt(20) = dCdt(20) - k_e2 * C_T7RNAP_M(j-1) * f_MM(j) * CTP_on(j); %CTP
        dCdt(21) = dCdt(21) - k_e2 * C_T7RNAP_M(j-1) * f_MM(j) * GTP_on(j); %GTP
        dCdt(23) = dCdt(23) + k_e2 * C_T7RNAP_M(j-1) * f_MM(j); %PPi
        dCdt(24) = dCdt(24) + k_e2 * C_T7RNAP_M(j-1) * f_MM(j); %H
    end

    dCdt(22) = 0; %Mg
end