function [obj,t_fine]=get_fine_field_t_all(obj,t,x_range,e_tx_out,w_0)
        %find the envelope peak.
        tau_fwhm_local= obj.tau*(sqrt(2*log(2)));
        dt=t(2)-t(1);
        dt_fine=2*pi/w_0/100;
        t_range_N=fix(obj.t_range/dt);
        %define the relevant time range for E_field ranging from the 
        %-fix(30e-15/dt) to 2t_N with respect to the peak position
        %the -fix(30e-15/dt) is the extra storage for spatial grating shift
        %later
        t_range_old=(0:2*t_range_N)+obj.t_start_px+fix(obj.ps_t/dt)...
            -fix(1.8*2*pi/w_0/dt)-0*fix(tau_fwhm_local/dt/2);
        %the fine time interpotation is symmetric with respect to 0
        t_interp_old=(-t_range_N:t_range_N)+fix(length(t)/2);
        t_fine=(t(fix(length(t)/2)-t_range_N):dt_fine:t(fix(length(t)/2)+t_range_N))';
        e_fine=zeros(length(t_fine),length(x_range));
        for kk=1:length(x_range)
        x_loc1=x_range(kk);
        et_slice_at_x=e_tx_out(t_range_old,x_loc1);
        e_fine(:,kk)= interp1(t(t_interp_old), et_slice_at_x,t_fine,'spline');
        end
        obj.e_fine=e_fine.*exp(1i.*w_0.*t_fine);

end
