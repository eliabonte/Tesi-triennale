function J = fun_lin(x)

k_tilde = [x.k1 x.k2 x.k3 x.k4]; 
assignin('base','k_tilde',k_tilde);

out = sim('controllo_leggelineare.slx');
output_xs2punti = getElement(out.yout, 'xs_2punti');
xs_2punti = output_xs2punti.Values.Data(:);


J = rms(xs_2punti);

end
