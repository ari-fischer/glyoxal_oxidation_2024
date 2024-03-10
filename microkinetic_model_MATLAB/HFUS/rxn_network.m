function rates = rxn_network(t,y,n_OH_nd,n_O_nd,ks_nd,ks_fit)
% input time, concentrations, OH rate, O rate, rate constants, rate
% constants fit

    global extents coeffs ind_OH ind_O2 ind_H2O C_gly0 ind_H3O l_ks ind_fit

    %assign values of parameters that are being regressed, modifying both
    %forward and reverse rates
    for i=1:length(ind_fit)
        k_scale = ks_fit(i)/ks_nd(ind_fit(i));
        ks_nd(ind_fit(i)) = ks_nd(ind_fit(i))*k_scale;
        ks_nd(ind_fit(i)+l_ks/2) = ks_nd(ind_fit(i)+l_ks/2)*k_scale;
    end
    
    %stack C vectors to get the C_matrix
    for i=1:length(ks_nd)
        c_M(:,i)=[y]; %each column is the concentration vector, each row is the same concentration repeated
    end
    %concentration matrix with rates to powers for reach step
    c_mij=c_M.^(coeffs');
    
    %rate vector with concentrations raised to their powers
    r_vector=ks_nd.*(prod(c_mij,1)');
    
    %get rates for species formation/consumption from reaction steps
    rates=extents'*r_vector;

    %assign OH to OH formation coefficient
    rates(ind_OH) = rates(ind_OH) + n_OH_nd;

    %assign O2 and H2O to zero 
    rates(ind_O2) = 0 ;
    rates(ind_H2O) = 0 ;
end