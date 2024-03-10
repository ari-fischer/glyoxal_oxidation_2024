function delRMSE = run_DEC(n_OH_fit, ks_nd, ks_fit,tspan,C0s_nd,ind_prod,ind_gly,ind_fit)
    %runs DRC analysis for small changed to rate constants of reaction steps
    % to determine which parameters contribute most to the error (Section
    % S2.9. Output the change in RMSE

    global ydata C_gly0 ts

    T_rxn = 42 %degrees C

    %get initial error by simulating rate product concentrations
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode23s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);

    %collect outputs
    outs_0 = get_error(ts,t,y,ind_prod,ydata,ind_gly);
    SSE_0 = sum(outs_0.^2);
    RMSE_0 = sqrt(mean(outs_0.^2));
    
    %initialize DRC calculations. Add one to the rate constants to reflect
    %OH formation rate
    DRC = ones([1,length(ks_nd)/2]);
    DRC = [DRC,1];
    delRMSE = []; % array of DRCs

    %use a constant perturbation with 0.01 change in free energy of TS
    pert = 0.01/8.3144E-3/(273+T_rxn);
    exp_pert = exp(-pert);
    
    for i = 1:length(DRC); %loop over DRC
        i %print step

        %assign perturbed values
        if i == length(DRC);
            n_OH_fit = n_OH_fit*exp(-pert);
        elseif any(ind_fit==i);
            ks_fit(find(ind_fit==i)) = ks_fit(find(ind_fit==i))*exp(-pert);
        else
            ks_nd(i) = ks_nd(i)*exp_pert;
            ks_nd(i+length(DRC)-1) = ks_nd(i+length(DRC)-1)*exp_pert;
        end

        %run simulation with updated values
        [t,y]=ode23s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);
        
        %recover initial value
        if i == length(DRC);
            n_OH_fit = n_OH_fit/exp(-pert);
        elseif any(ind_fit==i);
            ks_fit(find(ind_fit==i)) = ks_fit(find(ind_fit==i))/exp(-pert)
        else
            ks_nd(i) = ks_nd(i)/exp_pert;
            ks_nd(i+length(DRC)-1) = ks_nd(i+length(DRC)-1)/exp_pert;
        end
        
        %get error from outputs
        outs = get_error(ts,t,y,ind_prod,ydata,ind_gly);
    
        RMSE = sqrt(mean(outs.^2));
        delRMSE = [delRMSE,RMSE-RMSE_0];
    end

end