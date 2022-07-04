function [w,t,e_tx,e_fx]=E_field_analytic(obj,my_cons,electron,x)

            dx=x(2)-x(1);
            f=my_cons.f;
            A=my_cons.A;
            beta=my_cons.beta;
            w=my_cons.w_0+2*pi*((-125:0.08:125).*1e12)';
            delta_w=w-my_cons.w_0;
            k=my_cons.w_0/my_cons.c;
            df=(w(2)-w(1))/(2*pi);
            t = (-(length(w)-1)/(2*length(w)*df):1/(length(w)*df):(length(w)-1)/(2*length(w)*df))';
            dt=t(2)-t(1);
            s1=f;
            %f - grating size and the gap size
            d_from_focus=electron.g_h+electron.gap_h/2;
            s2=f-d_from_focus;
            M_11=A*(1-s2/f);
            M_12=(s1*(1-s2/f)+s2)/A;
            M_13=beta.*delta_w.*(s1*(1-s2/f)+s2);
            M_21=-A/f;
            M_22=(1-s1/f)/A;
            M_23=beta.*delta_w.*(1-s1/f);
            
            q_0=1i*pi*my_cons.sigma_0^2/my_cons.lambda_0;
            q_new=(M_11.*q_0+M_12)./(M_21.*q_0+M_22);
            
            e_fx=sqrt(1/(M_11+M_12/q_0)).*exp(-1i.*k.*(x-M_13).^2./q_new).*exp(-1i.*M_23.*k.*x)...
                .*exp(-my_cons.tau^2.*delta_w.^2/4).*exp(1i.*(my_cons.GDD.*delta_w.^2./2+my_cons.TOD.*delta_w.^3/6 ...
               +my_cons.FOD.*delta_w.^4/24 ));
            
            energy_sum=sum(sum(abs(e_fx).^2))*(0.5*my_cons.c.*my_cons.eps0)*df*dx;
            %normalized the amplitude in E_t to one
            e_fx=e_fx.*sqrt(obj.energy)/sqrt(energy_sum);
            e_tx  = fftshift(ifft(ifftshift(e_fx,1),[],1),1).*length(w)*df;

     
            
%             figure
%             imagesc(x*1e3,t./1e-12,abs(e_tx));
%             title(['|E(t,x)|, GDD=' num2str(my_cons.GDD*1e24) 'ps^2, TOD=' num2str(my_cons.TOD*1e36) 'ps^3'])
%             xlabel('x (mm)')
%             ylabel('t (ps)')
%             set(gca,'YDir','normal')
            
end
