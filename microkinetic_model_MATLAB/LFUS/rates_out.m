function out = rates_out(y,ks_nd,coeffs,extents,ks_fit,ind_fit) 
    %from concentrations of all species and their kinetic information,
    %calculate the rates for each elementary step and formation of species
    global l_ks
    for i=1:length(ind_fit)
        %use this to modify both forward and reverse rxn
        k_scale = ks_fit(i)/ks_nd(ind_fit(i));
        ks_nd(ind_fit(i)) = ks_nd(ind_fit(i))*k_scale;
        ks_nd(ind_fit(i)+l_ks/2) = ks_nd(ind_fit(i)+l_ks/2)*k_scale;
    end

    rates_matrix = []
    rates_species = []
    %follow calculation of rates from rxn_network function
    for j=1:length(y);
        for i=1:length(ks_nd);
            c_M(:,i)=[y(j,:)]; %each column is the concentration vector, each row is the same concentration repeated
        end
        %concentration matrix with rates to powers for reach step
        c_mij=c_M.^(coeffs');
        
        %rate vector with concentrations raised to their powers
        r_vector=ks_nd.*(prod(c_mij,1)');
        rates_matrix = [rates_matrix;r_vector'];
        rates_species = [rates_species;(extents'*r_vector)'];
    end

    out = [{rates_species},{rates_matrix}]
end