function out = rates_conv_plot(glys,rates_species,rates_matrix, n_OH_fit)
global C_gly0  ind_gly 
    %collect rates of important reactions to plot in figure 5 of main text

    conv = (C_gly0-glys(2:end))./C_gly0;
    y1 = rates_species*C_gly0*1e6; %uM/s
    y2 = rates_matrix*C_gly0*1e6;
    r_gly = y1(2:end,ind_gly);
    r_gly_OH = y2(2:end,6);
    r_glyox_n_OH = y2(2:end,10);
    r_glyox_OH = y2(2:end,8);
    r_ox_n_OH = y2(2:end,22);
    r_H2O2_n_OH = y2(2:end,5);

    r_OH = r_H2O2_n_OH.*0+n_OH_fit*C_gly0*1e6

    figure(20)
    hold on
    plot(conv,r_gly_OH+y2(2:end,28))
    plot(conv,r_gly_OH)
    plot(conv,r_glyox_n_OH)
    plot(conv,r_glyox_OH)
    plot(conv,r_ox_n_OH)
    plot(conv,r_H2O2_n_OH)
    plot(conv,r_OH)
    hold off
    out = [conv ,conv.*5,-r_gly,r_gly_OH,r_glyox_n_OH,...
        r_glyox_OH,r_ox_n_OH,r_H2O2_n_OH,r_OH]
    
    csvwrite("outputs/gly_rates.csv",out)
end