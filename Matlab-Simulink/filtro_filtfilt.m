function filtro = filtro_filtfilt(tipo_filtro,f_filter1,f_filter2,segnale)
filter = designfilt(tipo_filtro,'FilterOrder',4, ...
         'CutoffFrequency1',f_filter1,'CutoffFrequency2',f_filter2, ...
         'SampleRate',2);
%%filter=designfilt(tipo_filtro,'SampleRate',2,'HalfPowerFrequency',f_filter,'FilterOrder',4,...
        %%'DesignMethod','butter');
filtro=filtfilt(filter,segnale);    
end
