function Px_range = calcolorange_Pnn(uMR_values,f_max,p_molt)
%solo parte quadratica -> segnale da saturare (x^T * Px * x)
%con Px matrice simmetrica 4x4
%p11*x1^2 = fmax/2, 2*p12*x1*x2 = fmax/2
    
    [xspunto_max,xs_max,deltapunto_max,delta_max] = max_x(uMR_values);
    
    p11_max = (f_max/2)/(xspunto_max^2);
    p13_max = (f_max/2)/(xspunto_max*deltapunto_max);
    p14_max = (f_max/2)/(xspunto_max*delta_max);
    p33_max = (f_max/2)/(deltapunto_max^2);
    p34_max = (f_max/2)/(deltapunto_max*delta_max);
    p44_max = (f_max/2)/(delta_max^2);
    
    Px_range = p_molt*[p11_max p13_max p14_max p33_max p34_max p44_max];
end

