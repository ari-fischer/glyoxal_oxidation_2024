function out = run_ode(x) %get the t and y data, then the out is the difference
    global tspan C0s_nd ts ind_glyox_n ind_glyox C_gly0 ind_gly ind_ox ind_fit...
        ind_ox_n1 ind_ox_n2 ind_form ind_form_n ydata ks_nd_fit n_OH_0_fit ...
        ind_glyoxalate_dehyd ind_glyox_dehyd ind_gly_dehyd
    
    %assign coefficients with updated parameters
    n_OH_fit = n_OH_0_fit*x(1);
    ks_fit = [];
    %get the values to fit and push into the ODE solver
    for i=1:length(ind_fit)
        ks_fit = [ks_fit,ks_nd_fit(ind_fit(i))*x(i+1)];
    end
    ks_nd = ks_nd_fit;

    %simulate trends with rxn_network function
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode23s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
        ,tspan,C0s_nd);

    %% collect output concentrations from model
    % initialize reactant and product concentrations
    C_glyox_outs = [];
    C_ox_outs = [];
    C_form_outs = [];
    C_gly_outs = [];

    % collect the combined concentrations 
    for i=1:length(ts)
        t_sel = ts(i);
        sel_ind = sum(t<t_sel); %get the row for that data point
        
        %need to sum over reactants and products combining hydrates and
        %conjugate bases into same bins
        C_glyox_out = (y(sel_ind,ind_glyox)+y(sel_ind,ind_glyox_n) ...
            + y(sel_ind,ind_glyoxalate_dehyd) + y(sel_ind,ind_glyox_dehyd));
        
        C_ox_out = (y(sel_ind,ind_ox)+y(sel_ind,ind_ox_n1)+y(sel_ind,ind_ox_n2));

        C_form_out = (y(sel_ind,ind_form)+y(sel_ind,ind_form_n)); %find the column for glyox
        
        C_gly_out = y(sel_ind,ind_gly) + y(sel_ind,ind_gly_dehyd);

        C_glyox_outs = [C_glyox_outs;C_glyox_out];
        C_ox_outs = [C_ox_outs;C_ox_out];
        C_form_outs = [C_form_outs;C_form_out];
        C_gly_outs = [C_gly_outs;C_gly_out];
    end

    %calculate residuals renormalizing units
    out = ([C_ox_outs',C_glyox_outs',C_form_outs'].*C_gly0-ydata);
end