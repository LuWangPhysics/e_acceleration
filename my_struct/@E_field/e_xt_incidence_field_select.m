 function obj=e_xt_incidence_field_select(obj,t_fine)
                     dt=t_fine(2)-t_fine(1);
                      tau_fwhm_local= obj.tau*(sqrt(2*log(2)));
                      obj.t_start_p=find(real(obj.e_fine(:,2))==min(real(obj.e_fine(:,2))))...
                          +0*fix(tau_fwhm_local/dt/4);
                      dt=t_fine(2)-t_fine(1);
                      obj.t_start_p=obj.t_start_p+fix(obj.ps_t/dt);

                      obj.t=t_fine(obj.t_start_p:end);
                      obj.E_xt=obj.e_fine(obj.t_start_p:end,:);
 end
