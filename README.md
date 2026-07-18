%% ============================================================================
% Compact Model for Ferroelectric Negative Capacitance Transistor (MFIS-NCFET)
% Based on: Pahwa et al., IEEE TED, Vol. 64, No. 3, March 2017
% ============================================================================
% Clear workspace
clc; clear; close all;

%% Physical Constants
q = 1.602e-19;          % Elementary charge [C]
kb = 1.38e-23;          % Boltzmann constant [J/K]
eps0 = 8.854e-14;       % Vacuum permittivity [F/cm]
T = 300;                % Temperature [K]
Vt = kb * T / q;        % Thermal voltage [V]

%% Device Parameters (matching paper)
L = 1e-4;               % Channel length [cm] = 1 um
W = 1e-4;               % Channel width [cm] = 1 um
tox = 1e-7;             % Oxide thickness [cm] = 1 nm
eps_ox = 3.9;           % SiO2 relative permittivity
tfe = 5e-7;             % Ferroelectric thickness [cm] = 5 nm
Pr = 3e-6;              % Remnant polarization [C/cm^2] = 3 uC/cm^2
Ec = 1e6;               % Coercive field [V/cm] = 1 MV/cm
Na = 5e17;              % Body doping [cm^-3]
eps_s = 11.7;           % Si relative permittivity
phi_m = 4.1;            % Metal work function [eV]
phi_s = 4.05;           % Semiconductor work function [eV]
Eg = 1.12;              % Bandgap [eV]
mu = 500;               % Mobility [cm^2/Vs]

%% Derived Parameters
Cox = eps0 * eps_ox / tox;          % Oxide capacitance [F/cm^2]
eps_si = eps0 * eps_s;
gamma = sqrt(2 * q * eps_si * Na) / Cox;

ni = 1.5e10;                        % Intrinsic carrier concentration [cm^-3]
phi_F = Vt * log(Na / ni);          % Fermi potential [V]

VFB = phi_m - phi_s - Eg/2 - phi_F; % Flat band voltage [V]

% Landau coefficients: Vfe = a*QG + b*QG^3
% From Eq. (2) in paper:
% a = -3*sqrt(3)/2 * Ec*tfe/Pr
% b = 3*sqrt(3)/2 * Ec*tfe/Pr^3
a = -1.5 * sqrt(3) * Ec * tfe / Pr;
b = 1.5 * sqrt(3) * Ec * tfe / (Pr^3);
a_eff = a + 1.0 / Cox;

fprintf('Device Parameters:\n');
fprintf('  tfe = %.1f nm\n', tfe*1e7);
fprintf('  Pr = %.1f uC/cm^2\n', Pr*1e6);
fprintf('  Ec = %.1f MV/cm\n', Ec/1e6);
fprintf('  a = %.4e\n', a);
fprintf('  b = %.4e\n', b);
fprintf('  a_eff = %.4e\n', a_eff);
fprintf('  VFB = %.4f V\n', VFB);
fprintf('  Cox = %.4f uF/cm^2\n', Cox*1e6);
fprintf('  gamma = %.4f\n', gamma);
fprintf('  phi_F = %.4f V\n', phi_F);
fprintf('  2*phi_F = %.4f V\n', 2*phi_F);

%% ============================================================================
% Figure 2: Gate Charge, Surface Potential, Capacitance, Voltage Amplification
% ============================================================================
figure('Name', 'Figure 2: Gate Characteristics', 'Position', [100 100 1600 400]);

Pr_plot = 3e-6;
Ec_plot = 1e6;
VFB_plot = 0.0;
tfe_values = [0e-7, 1e-7, 3e-7, 5e-7];
tfe_labels = {'0nm', '1nm', '3nm', '5nm'};
% Use RGB triplets instead of character strings for colors
colors_rgb = [0 0 0;      % black for 0nm
              0 0.5 0;    % green for 1nm
              0.6 0.2 0;  % brown for 3nm
              0 0 1];     % blue for 5nm

VG_range = linspace(0, 2.0, 200);

% (a) Gate Charge Density QG vs VG
subplot(1, 4, 1);
for i = 1:length(tfe_values)
    [a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i] = ...
        get_params(tfe_values(i), Pr_plot, Ec_plot, VFB_plot, tox, eps_ox, eps_s, Na, T);
    
    QG_vals = zeros(size(VG_range));
    for j = 1:length(VG_range)
        QG_vals(j) = abs(calculate_QG(VG_range(j), 0.0, a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i));
    end
    plot(VG_range, QG_vals * 1e6, 'Color', colors_rgb(i,:), 'LineWidth', 2, 'DisplayName', tfe_labels{i});
    hold on;
end

% VC = 0.5V for tfe = 5nm
[a_5, b_5, a_eff_5, Cox_5, gamma_5, VFB_5, phi_F_5, Vt_5] = ...
    get_params(5e-7, Pr_plot, Ec_plot, VFB_plot, tox, eps_ox, eps_s, Na, T);
QG_vc = zeros(size(VG_range));
for j = 1:length(VG_range)
    QG_vc(j) = abs(calculate_QG(VG_range(j), 0.5, a_5, b_5, a_eff_5, Cox_5, gamma_5, VFB_5, phi_F_5, Vt_5));
end
plot(VG_range, QG_vc * 1e6, 'b.', 'MarkerSize', 4, 'DisplayName', '5nm, V_C=0.5V');

xlabel('V_G (V)', 'FontSize', 12);
ylabel('Q_G (\muC/cm^2)', 'FontSize', 12);
title('(a) Gate Charge Density', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
xlim([0 2.0]); ylim([0 3.5]);
grid on; box on;

% (b) Surface Potential psiS vs VG
subplot(1, 4, 2);
for i = 1:length(tfe_values)
    [a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i] = ...
        get_params(tfe_values(i), Pr_plot, Ec_plot, VFB_plot, tox, eps_ox, eps_s, Na, T);
    
    psiS_vals = zeros(size(VG_range));
    for j = 1:length(VG_range)
        psiS_vals(j) = calculate_psiS(VG_range(j), 0.0, a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i);
    end
    plot(VG_range, psiS_vals, 'Color', colors_rgb(i,:), 'LineWidth', 2, 'DisplayName', tfe_labels{i});
    hold on;
end

psiS_vc = zeros(size(VG_range));
for j = 1:length(VG_range)
    psiS_vc(j) = calculate_psiS(VG_range(j), 0.5, a_5, b_5, a_eff_5, Cox_5, gamma_5, VFB_5, phi_F_5, Vt_5);
end
plot(VG_range, psiS_vc, 'b.', 'MarkerSize', 4);

xlabel('V_G (V)', 'FontSize', 12);
ylabel('\psi_S (V)', 'FontSize', 12);
title('(b) Surface Potential', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
xlim([0 2.0]); ylim([0 1.6]);
grid on; box on;

% (c) Gate Capacitance CG vs VG
subplot(1, 4, 3);
for i = 1:length(tfe_values)
    [a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i] = ...
        get_params(tfe_values(i), Pr_plot, Ec_plot, VFB_plot, tox, eps_ox, eps_s, Na, T);
    
    CG_vals = zeros(size(VG_range));
    for j = 1:length(VG_range)
        CG_vals(j) = abs(calculate_CG(VG_range(j), 0.0, a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i));
    end
    plot(VG_range, CG_vals * 1e6, 'Color', colors_rgb(i,:), 'LineWidth', 2, 'DisplayName', tfe_labels{i});
    hold on;
end

CG_vc = zeros(size(VG_range));
for j = 1:length(VG_range)
    CG_vc(j) = abs(calculate_CG(VG_range(j), 0.5, a_5, b_5, a_eff_5, Cox_5, gamma_5, VFB_5, phi_F_5, Vt_5));
end
plot(VG_range, CG_vc * 1e6, 'b.', 'MarkerSize', 4);

xlabel('V_G (V)', 'FontSize', 12);
ylabel('C_G (\muF/cm^2)', 'FontSize', 12);
title('(c) Gate Capacitance', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0 2.0]); ylim([0 30]);
grid on; box on;

% (d) Voltage Amplification AV vs VG
subplot(1, 4, 4);
for i = 1:length(tfe_values)
    [a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i] = ...
        get_params(tfe_values(i), Pr_plot, Ec_plot, VFB_plot, tox, eps_ox, eps_s, Na, T);
    
    AV_vals = zeros(size(VG_range));
    for j = 1:length(VG_range)
        AV_vals(j) = calculate_AV(VG_range(j), 0.0, a_i, b_i, a_eff_i, Cox_i, gamma_i, VFB_i, phi_F_i, Vt_i);
    end
    plot(VG_range, AV_vals, 'Color', colors_rgb(i,:), 'LineWidth', 2, 'DisplayName', tfe_labels{i});
    hold on;
end

AV_vc = zeros(size(VG_range));
for j = 1:length(VG_range)
    AV_vc(j) = calculate_AV(VG_range(j), 0.5, a_5, b_5, a_eff_5, Cox_5, gamma_5, VFB_5, phi_F_5, Vt_5);
end
plot(VG_range, AV_vc, 'b.', 'MarkerSize', 4);

xlabel('V_G (V)', 'FontSize', 12);
ylabel('A_V = \partial\psi_S/\partialV_G', 'FontSize', 12);
title('(d) Voltage Amplification', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0 2.0]); ylim([0 3.0]);
grid on; box on;

sgtitle('Figure 2: Gate Characteristics of MFIS-NCFET', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 3: ID-VG and ID-VD Characteristics
% ============================================================================
figure('Name', 'Figure 3: Current Characteristics', 'Position', [50 50 1400 1000]);

[a_base, b_base, a_eff_base, Cox_base, gamma_base, VFB_base, phi_F_base, Vt_base] = ...
    get_params(0, Pr, Ec, -0.75, tox, eps_ox, eps_s, Na, T);
[a_nc, b_nc, a_eff_nc, Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc] = ...
    get_params(5e-7, Pr, Ec, -0.75, tox, eps_ox, eps_s, Na, T);

VG_range = linspace(0, 1.0, 150);
VD_range = linspace(0, 1.0, 150);

% (a) ID-VG characteristics
subplot(2, 2, 1);
for VD = [0.05, 1.0]
    ID_base = zeros(size(VG_range));
    for j = 1:length(VG_range)
        ID_base(j) = abs(calculate_ID(VG_range(j), VD, 0.0, a_base, b_base, a_eff_base, ...
            Cox_base, gamma_base, VFB_base, phi_F_base, Vt_nc, mu, W, L));
    end
    if VD == 0.05
        semilogy(VG_range, ID_base, '--', 'Color', [0.6 0.3 0], 'LineWidth', 2, 'DisplayName', 'Baseline V_D=50mV');
    else
        semilogy(VG_range, ID_base, '-.', 'Color', [0.6 0.3 0], 'LineWidth', 2, 'DisplayName', 'Baseline V_D=1V');
    end
    hold on;
end

for VD = [0.05, 1.0]
    ID_nc = zeros(size(VG_range));
    for j = 1:length(VG_range)
        ID_nc(j) = abs(calculate_ID(VG_range(j), VD, 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L));
    end
    if VD == 0.05
        semilogy(VG_range, ID_nc, 'r-', 'LineWidth', 2.5, 'DisplayName', 'NCFET V_D=50mV');
    else
        semilogy(VG_range, ID_nc, 'b-', 'LineWidth', 2.5, 'DisplayName', 'NCFET V_D=1V');
    end
end

% 60mV/dec reference
semilogy([0.15 0.35], [2e-4 2e-2], 'k--', 'LineWidth', 1.5, 'DisplayName', '60mV/dec');
text(0.32, 3e-2, '60mV/dec', 'FontSize', 9);

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(a) I_D - V_G Characteristics', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-6 1e3]);
grid on; box on;

% (b) ID-VD characteristics
subplot(2, 2, 2);
VG_values = [0.6, 0.8, 1.0];
colors_vg_rgb = [0 0.5 0; 1 0.6 0; 0 0 1];  % green, orange, blue

for i = 1:length(VG_values)
    VG = VG_values(i);
    ID_base = zeros(size(VD_range));
    ID_nc = zeros(size(VD_range));
    for j = 1:length(VD_range)
        ID_base(j) = abs(calculate_ID(VG, VD_range(j), 0.0, a_base, b_base, a_eff_base, ...
            Cox_base, gamma_base, VFB_base, phi_F_base, Vt_nc, mu, W, L));
        ID_nc(j) = abs(calculate_ID(VG, VD_range(j), 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L));
    end
    plot(VD_range, ID_base, '--', 'Color', colors_vg_rgb(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;
    plot(VD_range, ID_nc, '-', 'Color', colors_vg_rgb(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('V_G=%.1fV', VG));
end

xlabel('V_D (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(b) I_D - V_D Characteristics', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
xlim([0 1.0]); ylim([0 800]);
grid on; box on;

% (c) Transconductance gm vs VG
subplot(2, 2, 3);
for VD = [0.05, 1.0]
    gm_nc = zeros(size(VG_range));
    for j = 3:length(VG_range)-2
        dV = 1e-3;
        ID1 = calculate_ID(VG_range(j)-dV, VD, 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L);
        ID2 = calculate_ID(VG_range(j)+dV, VD, 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L);
        gm_nc(j) = abs((ID2 - ID1) / (2*dV) * 1e3); % mS
    end
    if VD == 0.05
        plot(VG_range, gm_nc, 'r-', 'LineWidth', 2.5, 'DisplayName', 'NCFET V_D=50mV');
    else
        plot(VG_range, gm_nc, 'b-', 'LineWidth', 2.5, 'DisplayName', 'NCFET V_D=1V');
    end
    hold on;
end

xlabel('V_G (V)', 'FontSize', 12);
ylabel('g_m (mS)', 'FontSize', 12);
title('(c) Transconductance', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
xlim([0 1.0]); ylim([0 1.2]);
grid on; box on;

% (d) Drain conductance gds vs VD
subplot(2, 2, 4);
for i = 1:length(VG_values)
    VG = VG_values(i);
    gds_nc = zeros(size(VD_range));
    for j = 3:length(VD_range)-2
        dV = 1e-3;
        ID1 = calculate_ID(VG, VD_range(j)-dV, 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L);
        ID2 = calculate_ID(VG, VD_range(j)+dV, 0.0, a_nc, b_nc, a_eff_nc, ...
            Cox_nc, gamma_nc, VFB_nc, phi_F_nc, Vt_nc, mu, W, L);
        gds_nc(j) = abs((ID2 - ID1) / (2*dV) * 1e3); % mS
    end
    semilogy(VD_range, gds_nc, '-', 'Color', colors_vg_rgb(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('V_G=%.1fV', VG));
    hold on;
end

xlabel('V_D (V)', 'FontSize', 12);
ylabel('g_{ds} (mS)', 'FontSize', 12);
title('(d) Drain Conductance', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-17 2.0]);
grid on; box on;

sgtitle('Figure 3: Current Characteristics of MFIS-NCFET', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 4: Ferroelectric Thickness Scaling
% ============================================================================
figure('Name', 'Figure 4: Thickness Scaling', 'Position', [100 100 1200 500]);

tfe_scales = [1e-7, 3e-7, 5e-7];
tfe_labels_s = {'1nm', '3nm', '5nm'};
colors_s_rgb = [0 0.5 0; 0.6 0.2 0; 0 0 1];  % green, brown, blue

% (a) ID-VG
subplot(1, 2, 1);
ID_base_s = zeros(size(VG_range));
for j = 1:length(VG_range)
    ID_base_s(j) = abs(calculate_ID(VG_range(j), 1.0, 0.0, a_base, b_base, a_eff_base, ...
        Cox_base, gamma_base, VFB_base, phi_F_base, Vt_nc, mu, W, L));
end
semilogy(VG_range, ID_base_s, 'k--', 'LineWidth', 2, 'DisplayName', 'Baseline FET');
hold on;

for i = 1:length(tfe_scales)
    [a_s, b_s, a_eff_s, Cox_s, gamma_s, VFB_s, phi_F_s, Vt_s] = ...
        get_params(tfe_scales(i), Pr, Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_s = zeros(size(VG_range));
    for j = 1:length(VG_range)
        ID_nc_s(j) = abs(calculate_ID(VG_range(j), 1.0, 0.0, a_s, b_s, a_eff_s, ...
            Cox_s, gamma_s, VFB_s, phi_F_s, Vt_s, mu, W, L));
    end
    semilogy(VG_range, ID_nc_s, 'Color', colors_s_rgb(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('t_{fe}=%s', tfe_labels_s{i}));
end

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(a) I_D - V_G (Thickness Scaling)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-4 1e3]);
grid on; box on;

% (b) ID-VD
subplot(1, 2, 2);
for i = 1:length(tfe_scales)
    [a_s, b_s, a_eff_s, Cox_s, gamma_s, VFB_s, phi_F_s, Vt_s] = ...
        get_params(tfe_scales(i), Pr, Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_s = zeros(size(VD_range));
    for j = 1:length(VD_range)
        ID_nc_s(j) = abs(calculate_ID(1.0, VD_range(j), 0.0, a_s, b_s, a_eff_s, ...
            Cox_s, gamma_s, VFB_s, phi_F_s, Vt_s, mu, W, L));
    end
    plot(VD_range, ID_nc_s, 'Color', colors_s_rgb(i,:), 'LineWidth', 2.5, 'DisplayName', sprintf('t_{fe}=%s', tfe_labels_s{i}));
    hold on;
end

xlabel('V_D (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(b) I_D - V_D (Thickness Scaling)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([0 800]);
grid on; box on;

sgtitle('Figure 4: Ferroelectric Thickness Scaling', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 5: Remnant Polarization Variation
% ============================================================================
figure('Name', 'Figure 5: Pr Variation', 'Position', [50 50 1600 400]);

Pr_values = [3e-6, 4e-6, 7e-6, 10e-6, 30e-6, 50e-6];
Pr_labels_v = {'3', '4', '7', '10', '30', '50'};
colors_pr = jet(length(Pr_values));

% (a) ID-VG
subplot(1, 4, 1);
for i = 1:length(Pr_values)
    [a_p, b_p, a_eff_p, Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p] = ...
        get_params(5e-7, Pr_values(i), Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_p = zeros(size(VG_range));
    for j = 1:length(VG_range)
        ID_nc_p(j) = abs(calculate_ID(VG_range(j), 1.0, 0.0, a_p, b_p, a_eff_p, ...
            Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p, mu, W, L));
    end
    semilogy(VG_range, ID_nc_p, 'Color', colors_pr(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('P_r=%s', Pr_labels_v{i}));
    hold on;
end
semilogy(VG_range, ID_base_s, 'k--', 'LineWidth', 2, 'DisplayName', 'Baseline');

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(a) I_D - V_G (P_r Variation)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 8);
xlim([0 1.0]); ylim([1e-6 1e3]);
grid on; box on;

% (b) ID-VD
subplot(1, 4, 2);
for i = 1:length(Pr_values)
    [a_p, b_p, a_eff_p, Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p] = ...
        get_params(5e-7, Pr_values(i), Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_p = zeros(size(VD_range));
    for j = 1:length(VD_range)
        ID_nc_p(j) = abs(calculate_ID(1.0, VD_range(j), 0.0, a_p, b_p, a_eff_p, ...
            Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p, mu, W, L));
    end
    plot(VD_range, ID_nc_p, 'Color', colors_pr(i,:), 'LineWidth', 2, ...
        'DisplayName', Pr_labels_v{i});
    hold on;
end

xlabel('V_D (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(b) I_D - V_D (P_r Variation)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 8);
xlim([0 1.0]); ylim([0 800]);
grid on; box on;

% (c) QGS-VG
subplot(1, 4, 3);
for i = 1:length(Pr_values)
    [a_p, b_p, a_eff_p, Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p] = ...
        get_params(5e-7, Pr_values(i), Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    QGS_vals = zeros(size(VG_range));
    for j = 1:length(VG_range)
        QGS_vals(j) = abs(calculate_QG(VG_range(j), 0.0, a_p, b_p, a_eff_p, ...
            Cox_p, gamma_p, VFB_p, phi_F_p, Vt_p));
    end
    plot(VG_range, QGS_vals * 1e6, 'Color', colors_pr(i,:), 'LineWidth', 2, ...
        'DisplayName', Pr_labels_v{i});
    hold on;
end

xlabel('V_G (V)', 'FontSize', 12);
ylabel('Q_{GS} (\muC/cm^2)', 'FontSize', 12);
title('(c) Q_{GS} - V_G', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 8);
xlim([0 1.0]); ylim([0 4.0]);
grid on; box on;

% (d) Loadline Analysis
subplot(1, 4, 4);
Vfe_range = linspace(-1.0, 1.0, 500);
for i = [1, 4, 6]  % Pr = 3, 10, 50
    [a_p, b_p, ~, ~, ~, ~, ~, ~] = ...
        get_params(5e-7, Pr_values(i), Ec, -0.75, tox, eps_ox, eps_s, Na, T);
    
    QG_landau = linspace(-15e-6, 15e-6, 500);
    Vfe_landau = a_p * QG_landau + b_p * QG_landau.^3;
    plot(Vfe_landau, QG_landau * 1e6, 'Color', colors_pr(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('P_r=%s', Pr_labels_v{i}));
    hold on;
end

% Loadlines
for VG_load = [0.2, 1.0]
    QG_load = zeros(size(Vfe_range));
    for j = 1:length(Vfe_range)
        psiS_approx = max(0, VG_load + 0.75 - Vfe_range(j));
        QG_load(j) = Cox * (VG_load + 0.75 - Vfe_range(j) - psiS_approx) * 1e6;
    end
    plot(Vfe_range, QG_load, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    text(Vfe_range(end-20), QG_load(end-20), sprintf('V_G=%.1fV', VG_load), 'FontSize', 9);
end

xlabel('V_{fe} (V)', 'FontSize', 12);
ylabel('Q_G (\muC/cm^2)', 'FontSize', 12);
title('(d) Loadline Analysis', 'FontSize', 13);
legend('Location', 'northwest', 'FontSize', 9);
xlim([-1.0 1.0]); ylim([-15 15]);
grid on; box on;
line([0 0], ylim, 'Color', 'k', 'LineWidth', 0.5);
line(xlim, [0 0], 'Color', 'k', 'LineWidth', 0.5);

sgtitle('Figure 5: Remnant Polarization Variation', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 6: Coercive Field Variation
% ============================================================================
figure('Name', 'Figure 6: Ec Variation', 'Position', [100 100 1200 500]);

Ec_values = [400e3, 600e3, 800e3, 1000e3];
Ec_labels_v = {'400', '600', '800', '1000'};
colors_ec_rgb = [0 0.5 0; 0 1 1; 1 0.6 0; 0 0 1];  % green, cyan, orange, blue

% (a) ID-VG
subplot(1, 2, 1);
for i = 1:length(Ec_values)
    [a_e, b_e, a_eff_e, Cox_e, gamma_e, VFB_e, phi_F_e, Vt_e] = ...
        get_params(5e-7, Pr, Ec_values(i), -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_e = zeros(size(VG_range));
    for j = 1:length(VG_range)
        ID_nc_e(j) = abs(calculate_ID(VG_range(j), 1.0, 0.0, a_e, b_e, a_eff_e, ...
            Cox_e, gamma_e, VFB_e, phi_F_e, Vt_e, mu, W, L));
    end
    semilogy(VG_range, ID_nc_e, 'Color', colors_ec_rgb(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('E_c=%skV/cm', Ec_labels_v{i}));
    hold on;
end
semilogy(VG_range, ID_base_s, 'k--', 'LineWidth', 2, 'DisplayName', 'Baseline FET');

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(a) I_D - V_G (E_c Variation)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-6 1e3]);
grid on; box on;

% (b) ID-VD
subplot(1, 2, 2);
for i = 1:length(Ec_values)
    [a_e, b_e, a_eff_e, Cox_e, gamma_e, VFB_e, phi_F_e, Vt_e] = ...
        get_params(5e-7, Pr, Ec_values(i), -0.75, tox, eps_ox, eps_s, Na, T);
    
    ID_nc_e = zeros(size(VD_range));
    for j = 1:length(VD_range)
        ID_nc_e(j) = abs(calculate_ID(1.0, VD_range(j), 0.0, a_e, b_e, a_eff_e, ...
            Cox_e, gamma_e, VFB_e, phi_F_e, Vt_e, mu, W, L));
    end
    plot(VD_range, ID_nc_e, 'Color', colors_ec_rgb(i,:), 'LineWidth', 2.5, ...
        'DisplayName', sprintf('E_c=%s', Ec_labels_v{i}));
    hold on;
end

xlabel('V_D (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title('(b) I_D - V_D (E_c Variation)', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([0 800]);
grid on; box on;

sgtitle('Figure 6: Coercive Field Variation', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 8: MFIS vs MFMIS Comparison
% ============================================================================
figure('Name', 'Figure 8: MFIS vs MFMIS', 'Position', [100 100 1200 500]);

[a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis] = ...
    get_params(5e-7, Pr, Ec, -0.75, tox, eps_ox, eps_s, Na, T);

% (a) Low VD = 50mV
subplot(1, 2, 1);
VD_low = 0.05;

ID_mfmis_low = zeros(size(VG_range));
ID_mfis_low = zeros(size(VG_range));
for j = 1:length(VG_range)
    % MFMIS approximation: constant Vint (VC=0)
    QG_mfmis = calculate_QG(VG_range(j), 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    psiS_mfmis = calculate_psiS(VG_range(j), 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    QI_mfmis = -QG_mfmis + gamma_mfis * Cox_mfis * sqrt(max(VG_range(j) - VFB_mfis - a_eff_mfis*QG_mfmis - b_mfis*QG_mfmis^3, 0));
    ID_mfmis_low(j) = abs(mu * (W/L) * QI_mfmis * VD_low);
    
    % MFIS: full model with spatial variation
    ID_mfis_low(j) = abs(calculate_ID(VG_range(j), VD_low, 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis, mu, W, L));
end

semilogy(VG_range, ID_mfmis_low, 'r-', 'LineWidth', 2.5, 'DisplayName', 'MFMIS');
hold on;
semilogy(VG_range, ID_mfis_low, 'b-', 'LineWidth', 2.5, 'DisplayName', 'MFIS');

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title(sprintf('(a) I_D - V_G at V_D = %.0fmV', VD_low*1000), 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-4 1e2]);
grid on; box on;

% (b) High VD = 1V
subplot(1, 2, 2);
VD_high = 1.0;

ID_mfmis_high = zeros(size(VG_range));
ID_mfis_high = zeros(size(VG_range));
for j = 1:length(VG_range)
    QG_mfmis = calculate_QG(VG_range(j), 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    psiS_mfmis = calculate_psiS(VG_range(j), 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    QI_mfmis = -QG_mfmis + gamma_mfis * Cox_mfis * sqrt(max(VG_range(j) - VFB_mfis - a_eff_mfis*QG_mfmis - b_mfis*QG_mfmis^3, 0));
    Vdsat = max(VG_range(j) - VFB_mfis - psiS_mfmis, 0.05);
    ID_mfmis_high(j) = abs(0.5 * mu * (W/L) * QI_mfmis * Vdsat);
    
    ID_mfis_high(j) = abs(calculate_ID(VG_range(j), VD_high, 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis, mu, W, L));
end

semilogy(VG_range, ID_mfmis_high, 'r-', 'LineWidth', 2.5, 'DisplayName', 'MFMIS');
hold on;
semilogy(VG_range, ID_mfis_high, 'b-', 'LineWidth', 2.5, 'DisplayName', 'MFIS');

xlabel('V_G (V)', 'FontSize', 12);
ylabel('I_D (\muA/\mum)', 'FontSize', 12);
title(sprintf('(b) I_D - V_G at V_D = %.0fV', VD_high), 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1e-4 1e2]);
grid on; box on;

sgtitle('Figure 8: MFIS vs MFMIS Comparison', 'FontSize', 14, 'FontWeight', 'bold');

%% ============================================================================
% Figure 9: Spatial Variation
% ============================================================================
figure('Name', 'Figure 9: Spatial Variation', 'Position', [100 100 1200 500]);

VG_fix = 1.0;
VD_values_sp = [0, 0.25, 0.5, 0.75, 1.0];
x_range = linspace(0, 1.0, 100);
colors_sp = parula(length(VD_values_sp));

% (a) QG(x) profiles
subplot(1, 2, 1);
for i = 1:length(VD_values_sp)
    VD = VD_values_sp(i);
    QG_S = calculate_QG(VG_fix, 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    QG_D = calculate_QG(VG_fix, VD, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    
    F_S = F_QG_func(QG_S, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis);
    F_D = F_QG_func(QG_D, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis);
    
    QG_profile = zeros(size(x_range));
    for j = 1:length(x_range)
        target_F = F_S - x_range(j) * (F_S - F_D);
        QG_profile(j) = fzero(@(QG) F_QG_func(QG, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis) - target_F, ...
            [min(QG_S, QG_D)-1e-6, max(QG_S, QG_D)+1e-6]) * 1e6;
    end
    
    plot(x_range, QG_profile, 'Color', colors_sp(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('V_D=%.2fV', VD));
    hold on;
end

xlabel('x/L', 'FontSize', 12);
ylabel('Q_G (\muC/cm^2)', 'FontSize', 12);
title('(a) Gate Charge Density Profile', 'FontSize', 13);
legend('Location', 'northeast', 'FontSize', 9);
xlim([0 1.0]); ylim([0 3.0]);
grid on; box on;

% (b) Vint(x) profiles
subplot(1, 2, 2);
for i = 1:length(VD_values_sp)
    VD = VD_values_sp(i);
    QG_S = calculate_QG(VG_fix, 0.0, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    QG_D = calculate_QG(VG_fix, VD, a_mfis, b_mfis, a_eff_mfis, ...
        Cox_mfis, gamma_mfis, VFB_mfis, phi_F_mfis, Vt_mfis);
    
    F_S = F_QG_func(QG_S, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis);
    F_D = F_QG_func(QG_D, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis);
    
    Vint_profile = zeros(size(x_range));
    for j = 1:length(x_range)
        target_F = F_S - x_range(j) * (F_S - F_D);
        QG_x = fzero(@(QG) F_QG_func(QG, VG_fix, a_mfis, b_mfis, a_eff_mfis, Cox_mfis, gamma_mfis, VFB_mfis) - target_F, ...
            [min(QG_S, QG_D)-1e-6, max(QG_S, QG_D)+1e-6]);
        Vfe_x = a_mfis * QG_x + b_mfis * QG_x^3;
        Vint_profile(j) = VG_fix - Vfe_x;
    end
    
    plot(x_range, Vint_profile, 'Color', colors_sp(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('V_D=%.2fV', VD));
    hold on;
end

line([0 1], [VG_fix VG_fix], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(0.85, VG_fix+0.05, 'V_G', 'FontSize', 9);
text(0.85, 2.25, 'NC', 'FontSize', 10, 'FontWeight', 'bold');
text(0.85, 1.95, 'PC', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('x/L', 'FontSize', 12);
ylabel('V_{int} (V)', 'FontSize', 12);
title('(b) Internal Voltage Profile', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 9);
xlim([0 1.0]); ylim([1.8 2.3]);
grid on; box on;

sgtitle('Figure 9: Spatial Variation Along Channel', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\nAll figures generated successfully!\n');