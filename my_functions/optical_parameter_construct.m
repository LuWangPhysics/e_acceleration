function my_cons=optical_parameter_construct(energy_in_x,beta_n)
my_cons=struct;
my_cons.c=3e8;
my_cons.eps0=8.854e-12;
%laser parameters 
my_cons.lambda_0=10e-6;
my_cons.w_0=my_cons.c*2*pi/my_cons.lambda_0;
my_cons.k_0=2*pi/my_cons.lambda_0;
my_cons.sigma_0=6e-3;


my_cons.tau_fwhm=100e-15;
my_cons.tau = my_cons.tau_fwhm/(sqrt(2*log(2))); 
my_cons.f_fwhm_total=0.44/my_cons.tau_fwhm;

%100fs
my_cons.GDD=4.5e-2*(1e-12)^2;
my_cons.TOD=0.94e-3*(1e-12)^3;

% my_cons.GDD=3.9e-2*(1e-12)^2;
% my_cons.TOD=1.1e-3*(1e-12)^3;



my_cons.FOD=0e-5*(1e-12)^4;

my_cons.f=100e-3;                                                          %focus lens
my_cons.A=1;                                                               % use littrow angle
%angular dispersion from grating
if abs(my_cons.GDD/2)>my_cons.tau^2/4
my_cons.beta=-beta_n*(my_cons.c/my_cons.w_0/my_cons.A/my_cons.sigma_0)*sqrt(abs(my_cons.GDD)-my_cons.tau^2/2);
else
    fprintf('minimal pulse duration cant be achieved');
end
end
