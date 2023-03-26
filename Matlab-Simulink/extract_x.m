function x = extract_x(uMR)
%faccio diverse simulazioni per diversi valori di uMR e per ogni
%simulazioni prendo il vettore delle variabili di stato 
%[xs_punto xs delta_punto delta]. Calcolo successivamente per ogni 
%variabile di stato il quantile al 95%. Essendo il segnale da saturare k*x e i
%limiti della saturazione [-fmax/2 fmax/2] i massimi valori utli sono
%kmax*xmax = fmax/2 --> invertendo la formula mi calcolo i valori del
%vettore k

assignin('base','val_uMR',uMR);
out = sim('QuarterCarModel_uMR.slx');

output_xspunto = getElement(out.yout, 'xs_punto');
xs_punto = output_xspunto.Values.Data(:);
qnt95_xspunto = quantile(abs(xs_punto),0.95); %calcolo il quantile al 95%

output_xs = getElement(out.yout, 'xs');
xs = output_xs.Values.Data(:);
qnt95_xs = quantile(abs(xs),0.95);

output_deltapunto = getElement(out.yout, 'delta_punto');
delta_punto = output_deltapunto.Values.Data(:);
qnt95_deltapunto = quantile(abs(delta_punto),0.95);

output_delta = getElement(out.yout, 'delta');
delta = output_delta.Values.Data(:);
qnt95_delta = quantile(abs(delta),0.95);

x = [qnt95_xspunto qnt95_xs qnt95_deltapunto qnt95_delta];
end

