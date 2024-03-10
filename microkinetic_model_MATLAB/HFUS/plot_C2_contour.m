function plot_C2_contour(n_OH_fit,ks_nd,ks_fit,C0s_nd,...
    coeffs,extents,ind_fit, time_solve)
global ind_glyox ind_glyox_n ind_glyox_dehyd ind_glyoxalate_dehyd...
    C_gly0 ind_ox ind_ox_n1 ind_ox_n2 ind_form ind_form_n ind_CO2...
     ind_bicarb ind_carbonic ind_gly ind_gly_dehyd l_ks ind_H3O

tspan = linspace(0,20000,1000);

get_pH = 1;
if get_pH == 1;
    %initialize data vectors
    gly_pH_M = [];
    glyox_pH_M = [];
    ox_pH_M = [];
    CO2_pH_M = [];
    form_pH_M = [];
    gly_conv_pH_M = [];
    pH_M = [];
    t_M = [];
    glyox_max = [];
    ox_max = [];
    C2s_max = [];

    gly_rates_M = [];
    glyox_rates_M = [];
    glyox_HT_M = [];
    oxal_HT_M = [];
    max_ratio = [];
    C_gly_max = [];
    t_max = [];
    glyoxylate_max = [];
    glyoxylic_acid_max = [];
    oxs_max_C2 = [];
    O2m_rxns_M = [];
    glyox_H2O2_rxns_M = [];

    pHs = linspace(0,5,6);
    %run constant pH simulations across range of values
    for i=1:length(pHs)
        pH_fix = pHs(i);
        %get concentration of H+
        C_H3O_fix = 10^(-pH_fix)/C_gly0;
        C0s_nd(ind_H3O) = C_H3O_fix;

        % run simulation
        options=odeset('RelTol',10^-6,'AbsTol',10^-6);
        [t,y]=ode23s(@(t,y) rxn_network_pH(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);

        %extract concentrations of products
        glyoxs = (y(:,ind_glyox)+y(:,ind_glyox_n)+y(:,ind_glyox_dehyd)+y(:,ind_glyoxalate_dehyd))*C_gly0;
        oxs = (y(:,ind_ox)+y(:,ind_ox_n1)+y(:,ind_ox_n2))*C_gly0;
        forms = (y(:,ind_form)+y(:,ind_form_n))*C_gly0;
        CO2s = (y(:,ind_CO2)+y(:,ind_bicarb)+y(:,ind_carbonic))*C_gly0;
        glys = (y(:,ind_gly)+y(:,ind_gly_dehyd))*C_gly0;

        glyoxylate = (y(:,ind_glyox_n)+y(:,ind_glyoxalate_dehyd))*C_gly0;
        glyoxylic = (y(:,ind_glyox)+y(:,ind_glyox_dehyd))*C_gly0;

        % get the rates of reaction steps and species out
        out = rates_out(y,ks_nd,coeffs,extents,ks_fit,ind_fit) ;
        rates_species=out{1};% rate of formation and consumption of each species
        rates_matrix = out{2};%rate of each elementary step

        %concentration of C2 acid products
        C2s = oxs+glyoxs;

        gly_conv = (C_gly0-glys)./C_gly0;

        %append matrices of outputs
        gly_pH_M = [gly_pH_M,glys];
        glyox_pH_M = [glyox_pH_M,glyoxs];
        ox_pH_M = [ox_pH_M,oxs];
        gly_conv_pH_M = [gly_conv_pH_M,gly_conv];
        CO2_pH_M = [CO2_pH_M,CO2s];
        form_pH_M = [form_pH_M,forms];
        pH_M = [pH_M,glys.*0+pH_fix];
        O2m_rxns_M = [O2m_rxns_M,rates_matrix(:,27)];
        glyox_H2O2_rxns_M = [glyox_H2O2_rxns_M,rates_matrix(:,19)];

        %maximum values of different specues
        t_M = [t_M,t];
        glyox_max=[glyox_max,max(glyoxs)];
        ox_max=[ox_max,max(oxs)];
        C2s_max = [C2s_max,max(C2s)];
        %concentrations at step that maximizes C2s
        C_gly_max = [C_gly_max,glys(find(C2s==max(C2s)))];
        max_ratio = [max_ratio,oxs(find(C2s==max(C2s)))./...
            glyoxs(find(C2s==max(C2s)))];

        glyoxylic_acid_max = [glyoxylic_acid_max,glyoxylic(find(C2s==max(C2s)))];
        glyoxylate_max = [  glyoxylate_max,glyoxylate(find(C2s==max(C2s)))];
        oxs_max_C2 = [oxs_max_C2,oxs(find(C2s==max(C2s)))]
        t_max = [t_max,t(find(C2s==max(C2s)))];

        %important reaction rates
        gly_rates_M = [gly_rates_M,rates_matrix(:,6)+rates_matrix(:,27)];
        glyox_rates_M = [glyox_rates_M,rates_species(:,ind_glyox_n)+rates_species(:,ind_glyox)]; 
        glyox_HT_M = [glyox_HT_M,rates_matrix(:,8)+rates_matrix(:,10)]; 

        oxal_HT_M = [oxal_HT_M,...
            rates_matrix(:,20)+rates_matrix(:,21)+rates_matrix(:,22)];
    end

end

% plot the maximum yields to C2 acids and other products at that reaction
% extent at different pH values, Figure 10
figure(10)
hold on
plot(pHs,C2s_max./C_gly0)
plot(pHs,C_gly_max./C_gly0)
plot(pHs,2*(C_gly0-C_gly_max-C2s_max)./C_gly0)
plot(pHs,C2s_max./(1+max_ratio)./C_gly0)%ox fac yield
plot(pHs,(C2s_max-C2s_max./(1+max_ratio))./C_gly0)%glyox fac yield
title('Maximum yield to C2 acids');
xlabel('pH');
ylabel('Fractional product yield');
legend('C2 acids','glyoxal','C1','oxalates','glyoxylates')
hold off

%plot the amount of OH formed at maximum C2 yield
figure(11)
plot(pHs,t_max*n_OH_fit.*C_gly0*1000)
title('OH consumed');
xlabel('pH');
ylabel('mM OH');

%export yields pHs, OH mM, fractional yield (of glyoxal reactant) to C2
%products, fractional yield to C1 products, oxalates yield, glyoxylates
%yield]
out1=[pHs',(t_max*n_OH_fit.*C_gly0*1000)',(C2s_max./C_gly0)',C_gly_max'./C_gly0,...
    (2*(C_gly0-C_gly_max-C2s_max)./C_gly0)',(C2s_max./(1+max_ratio)./C_gly0)'...
    (C2s_max-C2s_max./(1+max_ratio))'./C_gly0,t_max',...
    glyoxylic_acid_max'./C_gly0,glyoxylate_max'./C_gly0,...
    (glyoxylic_acid_max'+glyoxylate_max')./C_gly0,oxs_max_C2'./C_gly0];
csvwrite('outputs/C2_max_yield.csv',out1)


%generate contour plots Figures 9
X=t_M
Y=pH_M
Z1=glyox_pH_M
Z2=ox_pH_M
Z3=form_pH_M
Z4=CO2_pH_M

X_v = reshape(gly_conv_pH_M, [], 1); % Reshape X into a vector
Y_v = reshape(pH_M, [], 1); % Reshape Y into a vector
Z_v = reshape(glyox_pH_M, [], 1); % Reshape Z into a vector

%time in ks
ts_M = X./1000
%%
%plot the glyoxal/glyoxylic concetrations and rates
figure(12)
f = 10
c = 10
t_f = 30
hold on
tiledlayout(4,2)

%glyoxal concentrations
ax1 = nexttile;
[C, h1] = contourf(ts_M(2:end,:), Y(2:end,:), (C_gly0-gly_conv_pH_M(2:end,:).*C_gly0)*1E3,c,...
    'LineWidth', 1);
%title('glyoxal concentrations');
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
title({''})
xlim([0,t_f])
cb = colorbar;
h_round = round(h1.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');        % Add a color bar

ax2 = nexttile;

[C, h2] = contourf(ts_M(2:end,:), Y(2:end,:), gly_rates_M(2:end,:).*C_gly0.*1E6,c,'LineWidth', 1);
%title('glyoxal consumption rate');
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
title({''})
cb = colorbar;
h_round = round(h2.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');        % Add a color bar

ax3 = nexttile;
[C, h3] = contourf(ts_M(2:end,:), Y(2:end,:), glyox_pH_M(2:end,:)*1E3,c,'LineWidth', 1);
%title('glyoxylic acid concentration');
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
title({''})
cb = colorbar;
h_round = round(h3.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');        % Add a color bar

ax4 = nexttile;
[C, h4] = contourf(ts_M(2:end,:), Y(2:end,:), glyox_HT_M(2:end,:).*C_gly0.*1E6,c,'LineWidth', 1);
%title('glyoxylic acid HT rate');
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
title({''})
cb = colorbar;
h_round = round(h4.LevelList, 3) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');        % Add a color bar

ax5 = nexttile;
[C, h5] = contourf(ts_M(2:end,:), Y(2:end,:), ox_pH_M(2:end,:)*1E3,c,'LineWidth', 1);
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
title({''})
cb = colorbar;          % Add a color bar
h_round = round(h5.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');        % Add a color bar


ax6 = nexttile;
[C, h6] = contourf(ts_M(2:end,:), Y(2:end,:), oxal_HT_M(2:end,:).*C_gly0.*1E6,c,'LineWidth', 1);
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
title({''})
cb = colorbar;         % Add a color bar
h_round = round(h6.LevelList, 3) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');      % Add a color bar

ax7 = nexttile;
[C, h1] = contourf(ts_M(2:end,:), Y(2:end,:), ...
    (ox_pH_M(2:end,:)+glyox_pH_M(2:end,:)).*1E3,c,'LineWidth', 1);
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
cb = colorbar;         % Add a color bar
h_round = round(h1.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');      % Add a color bar

ax8 = nexttile;
[C, h2] = contourf(ts_M(2:end,:), Y(2:end,:), ...
    (form_pH_M(2:end,:)+CO2_pH_M(2:end,:)).*1E3,c,'LineWidth', 1);
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('pH','FontSize', f, 'FontName', 'Arial');
xlim([0,t_f])
cb = colorbar;         % Add a color bar
h_round = round(h2.LevelList, 2) %
% Set colorbar ticks to match contour levels
set(cb, 'Ticks', h_round,'FontSize', f, 'FontName', 'Arial');      % Add a color bar

set(gcf, 'Position', [100, 100, 550, 1000])
hold off
print('contour2','-dpng','-r300')

%save outputs as .csv files

Z1 = glyox_pH_M(2:end,:)*1E3
Z2 = glyox_HT_M(2:end,:).*C_gly0.*1E6
Z3 = ox_pH_M(2:end,:)*1E3
Z4 = oxal_HT_M(2:end,:).*C_gly0.*1E6
csvwrite('outputs/X.csv',X)
csvwrite('outputs/Y.csv',Y)
csvwrite('outputs/Z1.csv',Z1)
csvwrite('outputs/Z2.csv',Z2)
csvwrite('outputs/Z3.csv',Z3)
csvwrite('outputs/Z4.csv',Z4)

end