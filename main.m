clear all;
path(pathdef);
path_absolute = fileparts(which('main.m'));
addpath(path_absolute);
addpath([path_absolute '/my_functions/']);
addpath([path_absolute '/my_struct/']);

%------------------------------------
%laser parameters constants
%------------------------------------
save_flag=0;                                                               %1 save all the output figs
%scan N_shift from [-1,1] to find the maximum
N_shift=0;                                                              %N_shift*pi is the initial phases shift
beta_n=1;                                                                  %angular dispersion beta


name_head=['out' num2str(N_shift) 'shift_' num2str(beta_n) 'beta'];
sigma_y=0.2e-3;

energy=2.4*sigma_y;                                                        %[J]
energy_in_x=energy/sigma_y;                                                %[J/m]  
my_cons=optical_parameter_construct(energy_in_x,beta_n);
electron=electron_properties_construct(my_cons);

name_head=[num2str(my_cons.tau_fwhm.*1e15) 'fs'...
    num2str(fix(electron.energy_e/1e3/electron.e)) 'kev' name_head]; 


%------------------------------------
%construct electric field
%'tpf_stc_chirp'---stc scheme
%'direct'--- no tilt
%'tpf_standard_chirp'--- conventional pft
%------------------------------------
case_name= 'tpf_standard_chirp';
E=E_field;
E=E.initialize(energy_in_x,my_cons,electron,case_name,N_shift);

save_name=name_str(name_head,case_name,electron,E,my_cons,energy);
P=plot_set;
P=P.init(save_name,save_flag);

%------------------------------------
%construct input initial condition
%------------------------------------
NM=numerical_method;
[E,NM]=NM.init(my_cons,electron,E,case_name);
%for the perfect matched case
%[E,NM]=NM.init(my_cons,electron,E,'Euler_matched');


%when the time is longer than position mesh.
loop_break=@(v) E.x(end)-v*(E.t(2)-E.t(1));                               
for nn=2:length(E.t)
     if NM.z(nn-1)<loop_break(NM.v(nn-1))
     [E,NM]=NM.iteration(NM,nn,E);
     else
    
         break
     end

end


%-------------------------------------------------------------------
%calculate the instantaneous e field
%-------------------------------------------------------------------
 P.vz(E,NM,my_cons.c);
%-------------------------------------------------------------------
%calculate the kinetic energy gain 
%-------------------------------------------------------------------
 P.delta_kenetic(NM,my_cons,electron);

