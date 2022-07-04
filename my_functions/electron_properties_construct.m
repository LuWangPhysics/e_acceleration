function electron=electron_properties_construct(my_cons)

electron=struct;
electron.m=9.1e-31; 
electron.e=1.6e-19;  
electron.energy_e=20*1e3*electron.e;
electron.v_0=my_cons.c*sqrt(1-(1/(electron.energy_e/(electron.m*my_cons.c^2)+1))^2);
electron.tan_gamma=my_cons.c/electron.v_0;
electron.g_w=1*my_cons.lambda_0*(electron.v_0/my_cons.c);                          %spatial grating period
if electron.g_w<0.1e-6
    electron.g_w=1e-7;
end
electron.piller_w=electron.g_w/2;                                          %empty and the dielectric filled size is half of the period
electron.n_silica=3.4699;                                                  %spatial grating refractive index
electron.g_h=my_cons.lambda_0/2/(electron.n_silica-1);                     %spatial grating thickness
electron.gap_h=1e-6;

end