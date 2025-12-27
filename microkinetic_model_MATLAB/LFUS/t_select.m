function sel_inds = t_select(t,ts) 
    sel_inds = []
    for i=1:length(ts);
        t_sel = ts(i);
        sel_ind = sum(t<t_sel); %get the row for that data point
        sel_inds = [sel_inds,sel_ind]
    end
end