function [tan_tpf_angle,u,M,tau_fwhm_local,sigma_fwhm_total]=GDD_properties_at_focus(my_cons,beta,f,A)
q_0=1i*pi*my_cons.sigma_0^2/my_cons.lambda_0;
k=my_cons.w_0/my_cons.c;
q_new=-f^2/(q_0*A);
sigma_focus=sqrt(f^2/(abs(q_0)*A*k));
u=beta*f/sigma_focus^2;
M=my_cons.tau^2/4+u^2*sigma_focus^2;


tan_tpf_angle=atan(my_cons.c*u*my_cons.GDD/M)*180/pi;
tau_fwhm_local=2*sqrt(2*log(2)*(M+my_cons.GDD^2/4/M));
%f_fwhm_local=0.44/tau_fwhm_local;

sigma_fwhm_total=2*sqrt(log(2)*(sigma_focus^2/2+2*beta^2*f^2/my_cons.tau^2));
end