function out = get_error(ts,t,y,ind_prod,ydata,ind_gly)
    global C_gly0   
    ind_glyox = ind_prod(1);
    ind_glyox_n = ind_prod(2);
    ind_ox = ind_prod(3);
    ind_ox_n1 = ind_prod(4);
    ind_ox_n2 = ind_prod(5);
    ind_form = ind_prod(6);
    ind_form_n = ind_prod(7);
    ind_gly_dehyd = ind_prod(8);
    ind_glyox_dehyd = ind_prod(9);
    ind_glyoxalate_dehyd = ind_prod(10);
    
    %the follow code is constructed in the same way as in run_ODE

    C_glyox_outs = [];
    C_ox_outs = [];
    C_form_outs = [];
    C_gly_outs = [];

    for i=1:length(ts);
        t_sel = ts(i);
        sel_ind = sum(t<t_sel); %get the row for that data point
        
        %need to sum over
        C_glyox_out = (y(sel_ind,ind_glyox)+y(sel_ind,ind_glyox_n)...
            + y(sel_ind,ind_glyoxalate_dehyd) + y(sel_ind,ind_glyox_dehyd));
        
        C_ox_out = (y(sel_ind,ind_ox)+y(sel_ind,ind_ox_n1)+y(sel_ind,ind_ox_n2));

        C_form_out = (y(sel_ind,ind_form)+y(sel_ind,ind_form_n)); %find the column for glyox
        
        C_gly_out = y(sel_ind,ind_gly) + y(sel_ind,ind_gly_dehyd);

        C_glyox_outs = [C_glyox_outs;C_glyox_out];
        C_ox_outs = [C_ox_outs;C_ox_out];
        C_form_outs = [C_form_outs;C_form_out];
        C_gly_outs = [C_gly_outs;C_gly_out];
    end

    out = ([C_ox_outs',C_glyox_outs',C_form_outs'].*C_gly0-ydata);%;C_gly_outs]*C_gly0-ydata);

end