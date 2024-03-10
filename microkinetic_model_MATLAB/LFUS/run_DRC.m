function DRCs = run_DRC(n_OH_fit, ks_nd, ks_fit,tspan,C0s_nd,ind_prod,ind_gly,conv,ind_fit)
    global ydata C_gly0 ts
    % Now it is time to run DRC at 5% conversion or something to figure out
    % what the rate equation looks like for the formation of different
    % products. Use these to make the points about glyoxal consumption and
    % formic, oxalic, and glyoxylic acid formation rates.
    %initial rates
    ind_glyox = ind_prod(1);
    ind_glyox_n = ind_prod(2);
    ind_ox = ind_prod(3);
    ind_ox_n1 = ind_prod(4);
    ind_ox_n2 = ind_prod(5);
    ind_form = ind_prod(6);
    ind_form_n = ind_prod(7);
    ind_OH = ind_prod(8);

    
    options=odeset('RelTol',10^-6,'AbsTol',10^-6);
    [t,y]=ode15s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);

    glys = (y(:,ind_gly))*C_gly0;
    ind_sel = sum(((C_gly0-glys)./C_gly0)<conv)+1;
    %make new inits
    C0s_nd = y(ind_sel-1,:)
    time_solve = ts(end)*0.1
    t_nd = 1
    tspan = linspace(0,time_solve,200)./t_nd; 
    %calculate initial rates again
    
    [t,y]=ode15s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
            ,tspan,C0s_nd);

    rates0 = (y(2,:)-y(1,:))*C_gly0./(tspan(2)-tspan(1));
    %run simulation from 5% conversion
    
    r_glyox0 = rates0(ind_glyox)+rates0(ind_glyox_n);
    r_ox0 = rates0(ind_ox)+rates0(ind_ox_n1)+rates0(ind_ox_n2);
    r_form0 = rates0(ind_form)+rates0(ind_form_n);
    r_C2s0 = r_glyox0+r_ox0
    r_OH0 = rates0(ind_OH)

    DRC = ones([1,length(ks_nd)/2]);
    DRC = [DRC,1];
    Xs_gly = []; % array of DRCs
    Xs_glyox = [];
    Xs_ox = [];
    Xs_form = [];
    Xs_C2s = [];
    Xs_OH = [];


    
    pert = 0.001/8.314E-3/(315);
    exp_pert = exp(-pert);
    
    %need to add the c0s
    
    for i = 1:length(DRC); %loop over DRC
        i
        if i == length(DRC);
            n_OH_fit = n_OH_fit*exp(-pert);
        elseif any(ind_fit==i);
            ks_fit(find(ind_fit==i)) = ks_fit(find(ind_fit==i))*exp(-pert);
        else
            ks_nd(i) = ks_nd(i)*exp_pert;
            ks_nd(i+length(DRC)-1) = ks_nd(i+length(DRC)-1)*exp_pert;
        end
        [t,y]=ode15s(@(t,y) rxn_network(t,y,n_OH_fit,0,ks_nd,ks_fit)...
                ,tspan,C0s_nd);
        
        if i == length(DRC);
            n_OH_fit = n_OH_fit/exp(-pert);
        elseif any(ind_fit==i);
            ks_fit(find(ind_fit==i)) = ks_fit(find(ind_fit==i))/exp(-pert)
        else
            ks_nd(i) = ks_nd(i)/exp_pert;
            ks_nd(i+length(DRC)-1) = ks_nd(i+length(DRC)-1)/exp_pert;
        end
        % DRC at 5% conversion
        % find 5%?

        glys = (y(:,ind_gly))*C_gly0;
        rates = (y(2,:)-y(1,:))*C_gly0./(tspan(2)-tspan(1));
        %get rate of each species
    
        r_glyox = rates(ind_glyox)+rates(ind_glyox_n);
        r_ox = rates(ind_ox)+rates(ind_ox_n1)+rates(ind_ox_n2);
        r_form = rates(ind_form)+rates(ind_form_n);
        r_C2s = r_glyox + r_ox;
        r_OH =  rates(ind_OH)

        % are these at a particular conversion?
        X = (log(rates(ind_gly))-log(rates0(ind_gly)))/(-pert);
        Xs_gly = [Xs_gly,X]; 
        X = (log(r_glyox)-log(r_glyox0))/(-pert);
        Xs_glyox = [Xs_glyox,X];
        X = (log(r_ox)-log(r_ox0))/(-pert);
        Xs_ox = [Xs_ox,X];
        X = (log(r_form)-log(r_form0))/(-pert);
        Xs_form = [Xs_form,X];
        
        X = (log(r_C2s)-log(r_C2s0))/(-pert);
        Xs_C2s = [Xs_C2s,X];

        X = (log(r_OH)-log(r_OH0))/(-pert);
        Xs_OH = [Xs_OH,X];
    end
    %then add the nOH and nO
    
    
    %glyoxs = (y(:,ind_glyox)+y(:,ind_glyox_n))*C_gly0;
    %oxs = (y(:,ind_ox)+y(:,ind_ox_n1)+y(:,ind_ox_n2))*C_gly0;
    %forms = (y(:,ind_form)+y(:,ind_form_n))*C_gly0;
    
    sum(Xs_gly);
    sum(Xs_glyox);
    sum(Xs_ox);
    sum(Xs_form);
    sum(Xs_C2s);
    % table of all relevant DRC
    DRCs = [Xs_gly',Xs_glyox',Xs_ox',Xs_form',Xs_C2s',Xs_OH']

    


end