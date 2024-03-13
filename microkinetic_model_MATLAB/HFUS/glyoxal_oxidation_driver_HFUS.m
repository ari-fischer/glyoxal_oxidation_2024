%% 
clc
clear
clear all

%% Script to run regression and plotting for LFUS reactions at 42 C (315 K) 
% takes in experimental data, rate constants for reaction steps,
% stoichiometric coefficients, and reaction extents for the reaction
% network

%assign global variables to use in other functions

global tspan C0s_nd ts ind_glyox_n ind_glyox C_gly0 ind_gly ind_ox ind_fit... %run_ODE
        ind_ox_n1 ind_ox_n2 ind_form ind_form_n ydata ks_nd_fit n_OH_0_fit ... %run_ODE
        ind_glyoxalate_dehyd ind_glyox_dehyd ind_gly_dehyd ...%run_ODE
        extents coeffs ind_OH ind_O2 ind_H2O ind_H3O l_ks ind_OHn ... rxn_network
        ind_CO2 ind_bicarb ind_carbonic

    

%local path to folder containing experimental values "exp_data_HFUS.csv"
path_exp = "";
%path for input files generated with "get_export_ks_DFT_rxn_temp.ipynb"
%notebook
path_in = "";
exp_table = readtable(strcat(path_exp,'exp_data_HFUS.csv'));



%% Defining coefficients for the model

% initialize glyoxal, H+, OH- concentrations (moles/L)
C_gly0 = 5E-3;C_H3Op = 1E-7;C_OHn = 1E-7;C_H2O = 1;
t_nd = 1;


R_ig = 8.3144E-3; %kJ mol-1 K-1
kB = 1.380649E-23; % boltzmann's constant m2 kg s-2 K-1
h_pl = 6.62607015E-34; % planck's constant j s-1

% use henry's law for O2 concentration 
HL_O2 = 0.00099; % mol/kg/bar at 315 K
%O2 pressure and dissolved O2 concentration
P_O2 = 1; %bar
C_O2_0 = P_O2*HL_O2;

%reaction temperature (52 C for HFUS, 42 C for LFUS)
T = 52+273;

%%
% define kinetic and thermodynamic coefficients for reactions that are
% described with experimental, not DFT, values.

%for glyoxylate dehydration reactions
K_dehyd = 3.64E-3*exp(3.33E3/T); % equilibrium constant
k_dehyd = 2.04E9*exp(-7.45E3/T); %rate constant

%glyoxal dehydration kinetics/thermo at 298 K
T0_exp = 298;
% K = exp(-G/R/T)
K_dehyd_glyoxal = exp(T0_exp/T*log(350)); %equilibrium constant 
k_dehyd_glyoxal = 0.02; %rate constant

%oxalate reactions
k_oxalate_eT = 4.95E14*exp(-4.42E3/T);
k_Hoxalate_eT = 2.56E12*exp(-2.84E3/T);

%CO2 reactions
k_CO2_hydrate = 1.14E12*exp(-9.18E3/T);
K_CO2_hydrate = 2.62E3*exp(-1.76E3/T);

%OOH+OOH
k_2OOH = 8.4E5;

%H2O2_OH
k_H2O2_OH = 2.7E7;%

%OOH+OO-
k_OOH_O2n = 8.86E7;

% glyoxylate + H2O2 
k_glyoxn_H2O2 = 16.5;

% peroxyl disproportionation at 298 K
k_2OO = 4E8;

%proton transfer in water at 298 K
k_Hp_transfer = 1.4E11;

%hydrogen carbonate radical reaction with H2O2 reaction; 
k_carbonate_r_H2O2 = 4.3E5;

%% read network information
%read the reaction extents, coeffs, and rate constants generated from
%jupyter notebook 

%input extents of reaction from CSV
extents_table = readtable(strcat(path_in,'/rxn_extents.csv'));
%stoichiometric coefficients
coeffs_table = readtable(strcat(path_in,'/rxn_coeffs.csv'));
%rate constants
ks_table = readtable(strcat(path_in,'/rxn_ks.csv'));
%structures from SMILES
smiles_table = readtable(strcat(path_in,'/smiles.csv'));
species_cell = table2array(smiles_table(:,2));

%% define initial concentrations and match indices
C0s = ones(1,length(species_cell)).*0;

% get indices of important species by string comparisons
for i=1:length(species_cell);
    species_cell{i};
    if species_cell{i} == string('[OH]');
        ind_OH = i;
    elseif species_cell{i} == string('[OH][CH]([OH])[CH]([OH])[OH]');
        ind_gly = i;
        C0s(i) = C_gly0*(1-1/K_dehyd);
    elseif species_cell{i} == string('[O][O]')
        ind_O2 = i;
        C0s(i) = C_O2_0;
    elseif species_cell{i} == string('[OH2]')
        ind_H2O = i;
        C0s(i) = C_H2O;
    elseif species_cell{i} == string('[OH][C](=[O])[C](=[O])[OH]')
        ind_ox = i;
        %C0s(i) = C_gly0/10;
    elseif species_cell{i} == string('[OH][C](=[O])[C](=[O])[O-]')
        ind_ox_n1 = i;
    elseif species_cell{i} == string('[O-][C](=[O])[C](=[O])[O-]')
        ind_ox_n2 = i;
    elseif species_cell{i} == string('[OH][CH]([OH])[C](=[O])[OH]')
        ind_glyox = i;
    elseif species_cell{i} == string('[OH][CH]([OH])[C](=[O])[O-]')
        ind_glyox_n = i;
    elseif species_cell{i} == string('[CH]([OH])=[O]')
        ind_form = i;
    elseif species_cell{i} == string('[CH]([O-])=[O]')
        ind_form_n = i;
    elseif species_cell{i} == string('[O]')
        ind_O = i;
    elseif species_cell{i} == string('[O][OH]')
        ind_OOH = i;
    elseif species_cell{i} == string('[O][O-]')
        ind_OOn = i;
    elseif species_cell{i} == string('[OH][OH]')
        ind_H2O2 = i;
    elseif species_cell{i} == string('[OH3+]')
        ind_H3O = i;
        C0s(i) = C_H3Op;    
    elseif species_cell{i} == string('[OH-]');
        ind_OHn = i;
        C0s(i) = C_OHn;   
    elseif species_cell{i} == string('[OH][CH]([O][OH])[C]([OH])=[O]');
        ind_ROOH = i;
    elseif species_cell{i} == string('[O]=[C]=[O]')
        ind_CO2 = i;
    elseif species_cell{i} == string('[O][O][O]')
        ind_O3 = i;            %oxalic: [O]=[C]([OH])[C]([OH])=[O]
    elseif species_cell{i} == string('[OH][CH]([OH])[CH]=[O]')
        ind_gly_dehyd = i;      %glyoxylic: [O]=[CH][C]([OH])=[O]
        C0s(i) = C_gly0*(1/K_dehyd);
    elseif species_cell{i} == string('[O]=[CH][C](=[O])[OH]')
        ind_glyox_dehyd = i;      %glyoxylic: [O]=[CH][C]([OH])=[O]
    elseif species_cell{i} == string('[O]=[CH][C](=[O])[O-]')
        ind_glyoxalate_dehyd = i;      %glyoxylic: [O]=[CH][C]([OH])=[O]
    elseif species_cell{i} == string('[O-][C](=[O])[OH]')
        ind_bicarb = i;      %glyoxylic: [O]=[CH][C]([OH])=[O]
    elseif species_cell{i} == string('[O-][C](=[O])[O-]')
        ind_carbonic = i;      %glyoxylic: [O]=[CH][C]([OH])=[O]
    elseif species_cell{i} == string('[OH][CH]([OH])[C](=[O])[CH]=[O]')
        ind_propanal = i;
    elseif species_cell{i} == string('[OH][CH]([OH])[C]([O][O])([OH])[OH]')
        ind_peroxy1 = i;
    elseif species_cell{i} == string('[OH][C]([O][O])([OH])[C](=[O])[OH]')
        ind_peroxy2 = i;   
    elseif species_cell{i} == string('[OH][C]([O][O])([OH])[C](=[O])[O-]')
        ind_peroxy3 = i;    
    end

end

%index all of the important products
ind_prod = [ind_glyox,ind_glyox_n,ind_ox,...
    ind_ox_n1,ind_ox_n2,ind_form,ind_form_n,...
    ind_gly_dehyd, ind_glyox_dehyd, ind_glyoxalate_dehyd];



%access relevant data
%each column is a species, each row is an elementary step
%reaction extents for each species
extents = table2array(extents_table(:,3:end));
%stoichiometric coefficients
coeffs = table2array(coeffs_table(:,3:end));
orders = sum(coeffs,2); %nweeded for non-dimensionalization

%% add corrections to some of the rate constants

ks_0 = table2array(ks_table(:,end));
ks = ks_0;
%length of rate constants array
l_ks = length(ks);

ks(3) = k_OOH_O2n; % OOH O2- disprop 
ks(3+l_ks/2) = 0 %assume irreversible because reaction has large negative delta G

ks(5) = k_H2O2_OH; % H2O2 + OH reaction, 
ks(5+l_ks/2) = 0; %assume irreversible because reaction has large negative delta G

%kinetics for dehydration
%glyoxylic acid hydrate to aldehyde
ks(17) = k_dehyd;
ks(17+l_ks/2)= ks(17)*K_dehyd;
%glyxylate hydrate to aldehyde
ks(18) = k_dehyd;
ks(18+l_ks/2)= ks(18)*K_dehyd;

%glyoxal(2H2O) -> (H2O) 
ks(26) = k_dehyd_glyoxal;
ks(26+l_ks/2)= k_dehyd_glyoxal*K_dehyd_glyoxal;

%glyoxylic acid + H2O2 reaction
ks(19) = k_glyoxn_H2O2;
ks(19+l_ks/2) = 0; 

ks(21) = k_oxalate_eT;
ks(21+l_ks/2) = 0; 
ks(22) = k_Hoxalate_eT; 
ks(22+l_ks/2) = 0; 

%CO2 hydration
ks(23) = k_CO2_hydrate;
ks(23+l_ks/2) = k_CO2_hydrate/K_CO2_hydrate;

% Hydrogen carbonate radical H2O2 reaction
ks(25) = k_carbonate_r_H2O2 ;
ks(25+l_ks/2) = 0;

% lumped O2- addition to glyoxyl monohydrate and proton transfer 
ks(27) = 2E4*0.912544311;
ks(27+l_ks/2) = 0;

% O-transfer for peroxyl radicals from value reported for ethanol-derived
% species
ks(28) = k_2OO;

%non-dimensionalize all of the rate constants 
ks_nd = ks.*t_nd.*C_gly0.^(orders-1); % rate constants of the form s^-1 M^(n-1), where M is the reaction order


%% run the model

% get experimental data
ts = table2array(exp_table(:,1));
C_glyox = table2array(exp_table(:,5));
C_ox = table2array(exp_table(:,6));
C_form = table2array(exp_table(:,4));
C_gly_cons = table2array(exp_table(:,3));
C_gly = C_gly0 - C_glyox - C_ox
C_gly_cons_int = table2array(exp_table(:,7));

%compile into vectors
xdata = [ts',ts',ts'];
ydata = [C_ox',C_glyox',C_form'];

%time to solve the model for 
time_solve = ts(end)*1.5 ;
C0s_nd = C0s/C_gly0;

trials = length(C_ox);

% initialize ODE solver
tspan = linspace(0,time_solve,1000)./t_nd;

%volumetric OH formation rate; experiments-5.1000e-08
x0_exp = 5.1000e-08;

x0 = [2500e-6]/3600 ;
%non-dimensionalize OH formation
n_OH_0 = x0/C_gly0;

%initialize coefficients
beta = [0.690872609	4.711166251	0.344045585];

%specify indexes for rate coefficients to regress
ind_fit = [6 19];

%scale n_OH_0
n_OH_fit = n_OH_0*beta(1);

%run model to either fit beta, 
% 1 = regress beta; else model trends with initial values beta
run_ind = 1;

if run_ind == 1; %fit beta
    %specify initial values to send to regression function
    n_OH_0_fit = n_OH_0;
    ks_nd_fit = ks_nd;

    %bounds for parameter values
    lb = [0,0,0, 0];
    ub = [100,100,100, 100];

    %run regression with run_ode function
    options = optimoptions('lsqnonlin', 'DiffMinChange',.1,'Display','iter');
    [beta_fit,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(@run_ode,beta,lb,ub, options);
    ci = nlparci(beta_fit,residual,'jacobian',jacobian);
    
    %assign rate constants to regressed values
    n_OH_fit = n_OH_0*beta_fit(1);
    ks_fit = [];
    for i=1:length(ind_fit);
        ks_fit = [ks_fit,ks_nd_fit(ind_fit(i))*beta_fit(i+1)];
    end
    
    %simulate trends with regressed values using rxn_network function
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode23s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);

    %calculate regression statistics
    MAE = mean(abs(residual))*1000;

    RMSE = sqrt(mean(residual.^2))*1000;

    param_analy=1;
else %simulate trends with initial beta values
    
    %assign coefficients to values modified by beta
    n_OH_0_fit = n_OH_0;
    ks_nd_fit = ks_nd;
    n_OH_fit = n_OH_0*beta(1);
    ks_fit = []
    %get the values to fit and push into the ODE solver
    for i=1:length(ind_fit)
        ks_fit = [ks_fit,ks_nd_fit(ind_fit(i))*beta(i+1)];
    end
   
    %simulate trends rxn_network function
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode23s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);
    param_analy=0;
end

%%
%generate report for regressed data
%for table 1 in main text
if param_analy==1;
    param=[ci,beta_fit'];
    param(1,:) = param(1,:).*n_OH_0.*C_gly0;    
    param(2,:)=param(2,:).*ks(6)
    param(3,:)=param(3,:).*ks(19)

    mean_value_products = MAE/mean(ydata)/1000
end


%generate product concentrations from y matrix output
glyoxs = (y(:,ind_glyox)+y(:,ind_glyox_n)+y(:,ind_glyox_dehyd)+y(:,ind_glyoxalate_dehyd))*C_gly0;
oxs = (y(:,ind_ox)+y(:,ind_ox_n1)+y(:,ind_ox_n2))*C_gly0;
forms = (y(:,ind_form)+y(:,ind_form_n))*C_gly0;
CO2s = (y(:,ind_CO2)+y(:,ind_bicarb)+y(:,ind_carbonic))*C_gly0;
glys = (y(:,ind_gly)+y(:,ind_gly_dehyd))*C_gly0;
OHs = n_OH_0_fit*t*C_gly0*1000 %mM

ind_sel = t_select(t,ts)
glys(ind_sel)

%make conversion plot for C1 yields
model = [{glys},{glyoxs},{oxs},{forms},{CO2s}];
exp_cell = [{C_gly0-C_gly_cons_int},{C_glyox},{C_ox},{C_form}];

%% plot concentrations from experiments and simulations
figure(1)
hold on
f=12
plot(ts,C_glyox,'o',"MarkerEdgeColor","k",'MarkerSize',6,"MarkerFaceColor","#0072BD")
plot(ts,C_ox,'s',"MarkerEdgeColor","k",'MarkerSize',6,"MarkerFaceColor","#D95319")
plot(ts,C_form,'^',	"MarkerEdgeColor","k",'MarkerSize',6,"MarkerFaceColor","#77AC30")
%plot(ts,C_gly0-C_gly_cons,'d',	"MarkerEdgeColor","k",'MarkerSize',8,"MarkerFaceColor","k")
plot(ts,C_gly_cons,'d',	"MarkerEdgeColor","k",'MarkerSize',6,"MarkerFaceColor","k")
%lot(ts,C_glyox+C_ox+C_form,'x',"MarkerEdgeColor","k",'MarkerSize',8,"MarkerFaceColor","#0072BD")

plot(t,glyoxs,'--',"Color","#0072BD",'LineWidth',1)
plot(t,oxs,'--',"Color","#D95319",'LineWidth',1)
plot(t,forms,'--',"Color","#77AC30",'LineWidth',1)
plot(t,CO2s,'--',"Color","magenta",'LineWidth',1)
plot(t,C_gly0-glys,'--',"Color","k",'LineWidth',1)

%'ox acid-fit',
box on
xlabel("Time (s)",'FontSize',f, 'FontName', 'Arial')
ylabel("Yield (M)",'FontSize',f, 'FontName', 'Arial')
xlim([0 30000])
set(gcf, 'Position', [100, 100, 250, 280])
print('SI experiments and model nOH and kgly w O2n H2O2 SC_ODE23s','-dpng','-r300') %
hold off    
save('regress_full_model_ODE23s')

%% 
% run run_DEC function for error sensitivity discussed in Section S2.9
RMSE_calc = 1;
if RMSE_calc == 1;
    delRMSE = run_DEC(n_OH_fit, ks_nd, ks_fit,tspan,C0s_nd,ind_prod,ind_gly,ind_fit)
    %find(delRMSE==max(delRMSE));
    %max(delRMSE);
end


%outs for plotting figure 1 main text
[t(1:10:end)/1000,forms(1:10:end)*1000,glyoxs(1:10:end)*1000,...
    oxs(1:10:end)*1000,CO2s(1:10:end)*1000,C_gly0*1000-glys(1:10:end)*1000,...
    (C_gly0-glys(1:10:end))./C_gly0,glys(1:10:end)*1000]


%pH plot
figure(2)
plot(t/1000,-log10(y(:,ind_H3O)*C_gly0),'--',"Color","#0072BD",'LineWidth',1.5)
xlabel("Time (10^{3} s)",'FontSize',f, 'FontName', 'Arial')
ylabel("pH",'FontSize',f, 'FontName', 'Arial')

%pH and glyoxal/glyoxylic acid ratio for supporting information SI figure
%Section S2.8
figure(3)
hold on
f=12
subplot(1,2,1)
plot(t/1000,-log10(y(:,ind_H3O)*C_gly0),'--',"Color","#0072BD",'LineWidth',1.5)
xlabel("Time (10^{3} s)",'FontSize',f, 'FontName', 'Arial')
ylabel("pH",'FontSize',f, 'FontName', 'Arial')
%title('pH')
subplot(1,2,2)
plot(t(2:end)/1000,y(2:end,ind_glyox)./(y(2:end,ind_glyox)+y(2:end,ind_glyox_n)),'--',"Color","#0072BD",'LineWidth',1.5)
xlabel("Time (10^{3} s)",'FontSize',f, 'FontName', 'Arial')
ylabel("[glyox.H_{2}O]/[glyoxylates]",'FontSize',f, 'FontName', 'Arial')
ylim([0,1])
%title('[Glyox]/[Glyox-]')
hold off
set(gcf, 'Position', [100, 100, 400, 200])
print('SI time pH and glyoxylates','-dpng','-r300')

gly_cons = (C_gly0-glys)*1000;
figure(4)
hold on
f=12
subplot(1,2,1)
plot(gly_cons,-log10(y(:,ind_H3O)*C_gly0),'--',"Color","#0072BD",'LineWidth',1.5)
xlabel("Glyoxal reaction extent (mM)",'FontSize',f, 'FontName', 'Arial')
ylabel("pH",'FontSize',f, 'FontName', 'Arial')
%title('pH')
subplot(1,2,2)
plot(gly_cons(3:end),y(3:end,ind_glyox)./(y(3:end,ind_glyox)+y(3:end,ind_glyox_n)),'--',"Color","#0072BD",'LineWidth',1.5)
xlabel("Glyoxal reaction extent (mM)",'FontSize',f, 'FontName', 'Arial')
ylabel("[glyox.H_{2}O]/[glyoxylates]",'FontSize',f, 'FontName', 'Arial')
ylim([0,1])
%title('[Glyox]/[Glyox-]')
hold off
set(gcf, 'Position', [100, 100, 400, 200])
print('SI reaction extent pH and glyoxylates Section S2_8','-dpng','-r300')

%% figure S12 rates glyoxylates reactions 
% get the rates for elementary steps and for each species for section S2.8
out = rates_out(y,ks_nd,coeffs,extents,ks_fit,ind_fit) 
rates_species=out{1};
rates_matrix = out{2};

f=12
figure(5)
hold on
plot(t/1000,rates_matrix(:,8).*C_gly0.*1000000,'--','LineWidth',1.5)
plot(t/1000,rates_matrix(:,10).*C_gly0.*1000000,'--','LineWidth',1.5)
plot(t/1000,(rates_matrix(:,8)*ks(10)/ks(8)+rates_matrix(:,10)).*C_gly0.*1000000,'--','LineWidth',1.5)
xlabel("Time (10^{3} s)",'FontSize',f, 'FontName', 'Arial')
ylabel("Reaction rate ({\mu}M^{-1} s^{-1})",'FontSize',f, 'FontName', 'Arial')
%legend('glyoxylic','glyoxylate','without protonation')
box on
hold off
set(gcf, 'Position', [100, 100, 250, 280])
print('SI rates for glyoxylates reactions ','-dpng','-r300')

%average rate within the experimental time elapsed to report in text
r1 =rates_matrix(:,10).*C_gly0.*1000000;
r2 =(rates_matrix(:,8)*ks(10)/ks(8)+rates_matrix(:,10)).*C_gly0.*1000000;
mean(r2(1:ind_sel(end)))./mean(r1(1:ind_sel(end)));

%% table 1: predictions and benchmarks
%calculate the error for glyoxal consumed
ind_ts = [];
%find time indexs for experimental measurements
for i=1:length(ts)
    a = find(t>ts(i));
    ind_ts = [ind_ts,a(1)];
end

gly_cons = (C_gly0-glys)*1000;

mean_value = mean([C_glyox;C_ox;C_form])
if param_analy ==1
    %print for table 1
    %average glyoxal consumption rate from experiments
    A = mean(C_gly_cons([1,5,7,8])./ts([1,5,7,8]))
    %mean percentage error
    MAE/(mean_value*1000) % mean percentage error
    
    %gly % error
    gly_error = C_gly_cons.*1000-gly_cons(ind_ts);
    gly_error = gly_error(~isnan(gly_error));
    C_gly_clean = C_gly_cons(~isnan(C_gly_cons));
    mean(gly_error./(C_gly_clean.*1000))
   (C_gly_clean.*1000)- gly_error
end
%average gly consumption rate from model
mean(gly_cons(ind_ts)./t(ind_ts))

%% figure 5: rates of important reaction steps


%save the outputs for plotting as variable
out_fig4 = rates_conv_plot(glys,rates_species,rates_matrix,n_OH_fit);




%% fig 6: comparing formates and C2 acid acid yields

%experimental
[gly_cons(ind_ts),C_form,C_glyox+C_ox]
%model
[gly_cons,forms*1000,(glyoxs+oxs)*1000]




%% for fig 7 rate ratio
%get C2 rates
r_C2s_out = (rates_species(:,ind_glyox) + rates_species(:,ind_glyox_n)+ ...
    rates_species(:,ind_ox) + rates_species(:,ind_ox_n1) + rates_species(:,ind_ox_n2)).*C_gly0.*1000
%get forms rates
r_forms_out = (rates_species(:,ind_form) + rates_species(:,ind_form_n)).*C_gly0.*1000

%selectivity ratio
S_C2 = r_C2s_out./r_forms_out


%% fig 8: concentrations of important species 

%collect in matrix to plot in veusz
[t*n_OH_fit*C_gly0,(C_gly0-glys),y(:,ind_OH)*C_gly0.*1E9,y(:,ind_OOn).*C_gly0.*1000,...
    y(:,ind_glyox_n)*C_gly0,y(:,ind_H2O2)*C_gly0].*1000



%%
% peroxy radical concentrations at 10% conversion reported in SI
%ind_peroxy1, ind_peroxy2, ind_peroxy3
%a = find(gly_cons>C_gly0/10*1000);
%y(a(1),[ind_peroxy1, ind_peroxy2, ind_peroxy3])*C_gly0


%% Figure S9 reactions that consume OH
%list of reactions that involve OH
% 1, 4, 5, 6, 8, 10, 13, 15, 20, 21, 22, 24
%rates of elementary steps to export
% OH-OH; OH+OOH, OH+H2O2, gly+OH, glyox+OH, glyox- + OH, formic+OH,
% formate+OH, oxalic +OH, oxalate+OH, Hoxalate+OH, bicarbonate+OH

get_OHs = [1, 4, 5, 6, 8, 10, 13, 15, 20, 21, 22, 24]
rates_OH_out = rates_matrix(:,get_OHs)./n_OH_0_fit
[gly_cons,rates_OH_out]

%% acid base equilibrium Figure S3
% net rates
%net rate is forward - reverse
%approach to equilibrum is forward(1-rev/forward)
etas = rates_matrix(:,l_ks/2+1:end)./rates_matrix(:,1:l_ks/2);

figure(6)
tiledlayout(1,2)
ax1 = nexttile;
hold on
plot(gly_cons,etas(:,end-6))
plot(gly_cons,etas(:,end-5))
plot(gly_cons,etas(:,end-4))
plot(gly_cons,etas(:,end-3))
plot(gly_cons,etas(:,end-2))
plot(gly_cons,etas(:,end-1))
plot(gly_cons,etas(:,end))
f = 12
%legend('Glyox','Ox','Ox-','HCOOH','HCO3-','H2O','OOH','FontSize',f, 'FontName', 'Arial',...
%    'NumColumns',2)
ylim([.99,1.01])
xlim([0,0.1])
xlabel('Extent of glyoxal consumption ({\epsilon})','FontSize', f, 'FontName', 'Arial');
ylabel('Approach to equilibrium','FontSize', f, 'FontName', 'Arial');
ax = gca;
ax.FontSize = f;

hold off
ax2 = nexttile;
hold on
plot(t/1000,etas(:,end-6))
plot(t/1000,etas(:,end-5))
plot(t/1000,etas(:,end-4))
plot(t/1000,etas(:,end-3))
plot(t/1000,etas(:,end-2))
plot(t/1000,etas(:,end-1))
plot(t/1000,etas(:,end))
%legend('Glyox','Ox','Ox-','HCOOH','HCO3-','H2O','OOH','FontSize',f, 'FontName', 'Arial',...
%    'NumColumns',2)
ylim([.99,1.01])
xlabel('Time (10^{3} s)','FontSize', f, 'FontName', 'Arial');
ylabel('Approach to equilibrium','FontSize', f, 'FontName', 'Arial');
ax = gca;
ax.FontSize = f;
hold off
set(gcf, 'Position', [100, 100, 700, 300])
print('Approach to equilibrium','-dpng','-r300')

%% C2 contour plots; figure 9:
plot_C2_contour(n_OH_fit,ks_nd,ks_fit,C0s_nd,...
    coeffs,extents,ind_fit, time_solve)


%% DRC analysis to get the important steps Figure S4

%DEC analysis for each product
%use this to determine which parameters might need to be adjusted
convs = (C_gly0-glys)./C_gly0;
conv_span = linspace(0,.65,11);
out = rates_out(y,ks_nd,coeffs,extents,ks_fit,ind_fit) 
rates_species=out{1};
rates_matrix = out{2};
DRCs_formates = [];
DRCs_gly = [];
DRCs_glyox = [];
DRCs_ox = [];
DRCs_C2 = [];
DRCs_OH = [];

%DRC for epsilon LH
conv_span = linspace(0,0.65/5,11);

for i=1:length(conv_span)-1
    ind_sel = sum(((C_gly0-glys)./C_gly0)<conv_span(i+1));
    rates_conv = rates_matrix(ind_sel,:)
    DRCs = run_DRC(n_OH_fit, ks_nd, ks_fit,tspan,C0s_nd,[ind_prod,ind_OH],ind_gly,conv_span(i+1),ind_fit)
    
    %DRCS for each group of products
    DRCs_formates = [DRCs_formates,DRCs(:,4)];
    DRCs_gly = [DRCs_gly,DRCs(:,1)];
    DRCs_glyox = [DRCs_glyox,DRCs(:,2)];
    DRCs_ox = [DRCs_ox,DRCs(:,3)];
    DRCs_C2 = [DRCs_C2,DRCs(:,5)];
    DRCs_OH = [DRCs_OH,DRCs(:,6)];
end

%combining to copy into veusz
[conv_span(2:end)'*C_gly0*1000,DRCs_C2']
[conv_span(2:end)'*C_gly0*1000,DRCs_formates']
[conv_span(2:end)'*C_gly0*1000,DRCs_C2'-DRCs_formates']



