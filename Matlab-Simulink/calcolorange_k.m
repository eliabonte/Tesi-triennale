function k_vec = calcolorange_k(uMR_values,f_max,k_molt)
%Essendo il segnale da saturare k*x e i
%limiti della saturazione [-fmax/2 fmax/2] i massimi valori utli sono
%kmax*xmax = fmax/2 --> invertendo la formula mi calcolo i valori del
%vettore k

    [xspunto_max,xs_max,deltapunto_max,delta_max] = max_x(uMR_values);
    
    k1_max = (f_max/2)/(xspunto_max);
    k2_max = (f_max/2)/(xs_max);
    k3_max = (f_max/2)/(deltapunto_max);
    k4_max = (f_max/2)/(delta_max);
    
    k_vec = k_molt*[k1_max k3_max k2_max k4_max];
end

