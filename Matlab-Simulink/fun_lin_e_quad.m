function J = fun_lin_e_quad(x)

k_tilde = [x.k1 0 x.k3 x.k4];
Px = [0 0 x.p13  x.p14; 0 0 0 0; x.p13 0 0  x.p34;  x.p14 0  x.p34 0];

assignin('base','k_tilde',k_tilde);
assignin('base','Px',Px);

out_linquad = sim('controllo_leggelin_e_quad.slx');
output_xs2punti_linquad = getElement(out_linquad.yout, 'xs_2punti');
xs_2punti_linquad = output_xs2punti_linquad.Values.Data(:);

J = rms(xs_2punti_linquad);

end
