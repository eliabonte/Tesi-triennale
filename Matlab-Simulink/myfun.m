function J_sim = myfun(x)

xval = x.uMR; 
assignin('base','val_uMR',xval);

out = sim('QuarterCarModel_uMR.slx');
output_xs2punti = getElement(out.yout, 'xs_2punti');
xs_2punti = output_xs2punti.Values.Data(:);
    
J_sim = rms(xs_2punti);

end