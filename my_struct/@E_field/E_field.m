classdef E_field
    properties
        E0;
        tau;
        w_0;
        c;
        piller_w;
        v_0;
        t;
        x;
        t_plot;
        e;
        m;
        beam_fwhm,z_center,n_gaussian;
        v_analytic;
        tan_tpf_angle;
        energy;
        fluence;
        E_xt;
        field;
        E_z_electron;
        e_fine;
        e_gross;
        N_shift;
       %tpf_spatial_chirp_t_start;
      
       t_gross;
       t_start_p;
       t_start_px;   %the location of position x find t maximum
       t_range;
       ps_t;
       E_peak_t;
       d_lambda;
    end
methods
                [w,t,e_tx,e_fx]=E_field_analytic(obj,my_cons,electron,x_arr);
                [w,t,e_tx,e_fx]=E_field_tpf_standard(obj,my_cons,electron,x_arr);

                [obj,t_fine]=get_fine_field_t_all(obj,t,t_lim,x_range,e_tx_out,w_0);
                [obj,x_start,x_range]=e_x_select(obj,x_arr,e_xt_out);   
                obj=e_xt_incidence_field_select(obj,t_fine);
                obj=E_peak_t_each_x(obj,z,nn,NM);
                
                function [obj]=initialize(obj,energy,my_cons,electron,case_name,N_shift)
                    obj.N_shift=N_shift;
                    obj.m=electron.m;
                    obj.e=electron.e;
                    obj.piller_w=electron.piller_w;
                  
                    obj.energy=energy;
                    obj.w_0=my_cons.w_0;
                    obj.c=my_cons.c;
                    obj.v_0=electron.v_0;

                    %obj.beam_size=my_cons.sigma_0;
                    obj.n_gaussian=1;
                    if energy >2
                        obj.t_range=0.9e-12;
                    else
                    obj.t_range=1.4e-12;
                    end
                    obj.t=-obj.t_range:1/my_cons.w_0/25:obj.t_range;   %[s]
                 
                    %select the laser field
                    obj=obj.field_select(case_name,my_cons,electron);
                
                 
                    obj.E_z_electron=zeros(size(obj.t));
                    obj.E_peak_t=zeros(size(obj.t));
                 
                end
                
                %---------------------------------------------------------
                %define the type of electric field based on input
                %---------------------------------------------------------

                function obj=field_select(obj,case_name,my_cons,electron)
                    
                         %------------------------------
                         %get the x array and beam size
                         %------------------------------
                         [obj.tan_tpf_angle,~,~,~,obj.beam_fwhm]=GDD_properties_at_focus(my_cons,my_cons.beta,my_cons.f,my_cons.A);
                     
                         obj.z_center=0.45*obj.beam_fwhm;
                         x_arr=-obj.beam_fwhm*2.8:my_cons.lambda_0/300:obj.beam_fwhm*2.8;
                         %------------------------------
                         %get the equivalent tau_fwhm
                         %------------------------------
                         [~,t_arr,~,e_fx]=E_field_analytic(obj,my_cons,electron,x_arr); 
                         dt=t_arr(2)-t_arr(1);
                         e_xt_out=fftshift(ifft(ifftshift(e_fx,1),[],1),1)./dt;
                         [tau_fwhm_local]=funfwhm(t_arr,abs(e_xt_out(:,fix(length(x_arr)/2))).^2);
                          obj.tau=tau_fwhm_local/(sqrt(2*log(2))); 
                         
                         %phase shifter at t
                         %for 2.4 high energy the ps_t=0 standard, 0.5 chirp
                         %for low energy ps_t=0.75 for standard, 0 chirp
                         obj.ps_t=obj.N_shift*(2*pi/my_cons.w_0);
                         sigma=obj.beam_fwhm/sqrt(2*log(2));

                         switch case_name

  
                            case 'direct'
                                 [w,t_arr,e_fx]=E_field_direct(obj,my_cons,electron,x_arr);
                 
                                 dt=t_arr(2)-t_arr(1);
                                  e_xt_out=fftshift(ifft(ifftshift(e_fx,1),[],1),1)./dt;
                                  obj.E0=max(max(abs(e_xt_out)));
                                  %shift t postion
                                  [obj,x_start,x_range]=e_x_select(obj,x_arr,e_xt_out);   
                                  
                                  [obj,t_fine]=get_fine_field_t_all(obj,t_arr,x_range,e_xt_out(:,x_start:end),my_cons.w_0);
                                   obj=e_xt_incidence_field_select(obj,t_fine);
                                   obj.field=@numerical;

                             

                             case 'tpf_standard_chirp'  
                                 
                                  [w,t_arr,e_fx,obj.tan_tpf_angle]=E_field_tpf_standard(obj,my_cons,electron,x_arr);
                                   dt=t_arr(2)-t_arr(1);
                                   e_xt_out=fftshift(ifft(ifftshift(e_fx,1),[],1),1)./dt;
                                   obj.E0=max(max(abs(e_xt_out)));
                                  [obj,x_start,x_range]=e_x_select(obj,x_arr,e_xt_out);   
                                  
                                  [obj,t_fine]=get_fine_field_t_all(obj,t_arr,x_range,e_xt_out(:,x_start:end),my_cons.w_0);
                                   obj=e_xt_incidence_field_select(obj,t_fine);
                                   obj.field=@numerical;
                                   obj.e_gross=e_xt_out(:,x_start:end);
                             case 'tpf_stc_chirp'
                                  [w,t_arr,~,e_fx]=E_field_analytic(obj,my_cons,electron,x_arr); 
                                   dt=t_arr(2)-t_arr(1);
                                   %(t,x)
                                   e_xt_out=fftshift(ifft(ifftshift(e_fx,1),[],1),1)./dt;
                                   obj.fluence=max(my_cons.c*my_cons.eps0*0.5*sum(abs(e_xt_out).^2,1)*dt);
                                   obj.E0=max(max(abs(e_xt_out)));
                                   [obj,x_start,x_range]=e_x_select(obj,x_arr,e_xt_out); 
                                  
                                   [obj,t_fine]=get_fine_field_t_all(obj,t_arr,x_range,e_xt_out(:,x_start:end),my_cons.w_0);
                                   %get the negative maximum field strength
                                   %for input inital
                                   obj=e_xt_incidence_field_select(obj,t_fine);
                                   obj.field=@numerical;
                                   obj.d_lambda=lambda_x(obj,my_cons);
                                   obj.e_gross=e_xt_out(:,x_start:end);
                                   
                                      
                                 
                        end
  
                end
                 
                function d_lambda=lambda_x(obj,my_cons)
                    q_0=pi*my_cons.sigma_0^2/my_cons.lambda_0;
                    k=my_cons.w_0/my_cons.c;
                    q2=abs(-my_cons.f^2/q_0);
                    x_eff=obj.x-obj.z_center;

                    f_x=(my_cons.w_0+k*my_cons.beta*my_cons.f.*x_eff./(my_cons.tau.^2/4*q2+my_cons.beta^2*my_cons.f^2*k))./(2*pi);

                    d_lambda=my_cons.c./f_x;
                end


        
    end
end 