function J = fun_quad(x)

%%k_tilde = [x.k1 x.k2 x.k3 x.k4];
Px = [x.p11 0 x.p13  x.p14; 0 0 0 0; x.p13 0 x.p33  x.p34;  x.p14 0  x.p34 x.p44];

%%assignin('base','k_tilde',k_tilde);
assignin('base','Px',Px);

out_quad = sim('controllo_leggequadratica.slx');
output_xs2punti_quad = getElement(out_quad.yout, 'xs_2punti');
xs_2punti_quad = output_xs2punti_quad.Values.Data(:);

J = rms(xs_2punti_quad);

end
