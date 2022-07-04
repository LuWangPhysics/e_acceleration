function [E,obj]=Euler_chirp(obj,nn,E)
    

      obj=obj.Euler_z(nn);
                 
      if (sum(abs(obj.z(nn)-obj.z(1)))>obj.g_w_sum)||(nn==2)
       N_margin=3;
      [x_range,t_period,obj]=g_w_period_construct(obj,E,nn,N_margin);
      if (obj.N_x+x_range)>length(E.x)
%        x=0:dx:(length(E.x)-obj.N_x-1)*dx;
%        x_range=length(x);
        return;
      end
      
      dt=E.t(2)-E.t(1);
      t_range=fix(t_period/dt)+1;
      
      obj.N_x=find_x_N(E.x,obj.z(nn));
      
      %N_x=fix(obj.z(nn-1)/dx);
      if (obj.N_x+x_range)>length(E.x)
        return;
      end
     
      for kk=1:x_range
        x_i=obj.N_x+kk-1;
        t_start=E.t_start_p-t_range(kk);
        t_end=t_start+length(E.t)-1;
        E.E_xt(1:end,x_i)=E.e_fine(t_start:t_end,x_i);
      end
      obj.g_w_sum=obj.g_w_sum+obj.g_w_period(obj.g_w_n);
     
      obj.x_range=x_range-N_margin;
      obj.N_x=obj.N_x+obj.x_range;
      obj.test(obj.g_w_n)=obj.N_x;
      obj.g_w_n=obj.g_w_n+1;
      end
            

                 %calculate force using B and E at mesh 1,2,3,4   
                  %add evanancent field
                 k_z_decay=sqrt((2*pi/obj.g_w_period(obj.g_w_n-1))^2-(2*pi/obj.lambda_0).^2);
                 A_evan=exp(-k_z_decay*obj.gap_h/2);
                 obj.field_inst=A_evan.*E.field(E,obj.z(nn),E.t(nn));
                 %obj.field_inst=E.field(E,obj.z(nn),E.t(nn));
                 obj=obj.Euler_v(E,nn,real(obj.field_inst));
                 E.E_z_electron(nn)=obj.field_inst;
end