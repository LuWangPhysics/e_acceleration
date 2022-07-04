function [w,t,e_fx,tpf]=E_field_tpf_standard(obj,my_cons,electron,x)
           
            k=my_cons.w_0/my_cons.c;
          
            tpf=electron.tan_gamma;
            beta=tpf/my_cons.c/k;
            w=my_cons.w_0+2*pi*((-125:0.08:125).*1e12)';
            delta_w=w-my_cons.w_0;

            df=(w(2)-w(1))/(2*pi);
            t = (-(length(w)-1)/(2*length(w)*df):1/(length(w)*df):(length(w)-1)/(2*length(w)*df))';
            dt=t(2)-t(1);
            dx=x(2)-x(1);
            sigma=obj.beam_fwhm/sqrt(2*log(2));
       
          
            
            e_fx=exp(-x.^2./sigma^2).*exp(-1i.*beta.*delta_w.*k.*x).*exp(-obj.tau^2.*delta_w.^2/4);

            energy_sum=sum(sum(abs(e_fx).^2))*(0.5*my_cons.c.*my_cons.eps0)*df*dx;
            %normalized the amplitude in E_t to one
            e_fx=e_fx.*sqrt(obj.energy)/sqrt(energy_sum);
         
%             figure
%             imagesc(x*1e3,t./1e-12,abs(e_tx));
%             title(['|E(t,x)|, GDD=' num2str(my_cons.GDD*1e24) 'ps^2, TOD=' num2str(my_cons.TOD*1e36) 'ps^3'])
%             xlabel('x (mm)')
%             ylabel('t (ps)')
%             set(gca,'YDir','normal')
            
end
