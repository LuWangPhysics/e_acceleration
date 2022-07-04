classdef numerical_method
    properties
        dt,dz;
        p_array;
        v,z;
        field_inst;

        iteration;
        g_h;
        n_silica;
        c;
        g_w_sum;
        lambda_0;
        g_w_period;
        N_x;
        x_range;
        g_w_n;
        gap_h;
        E_ave;
        W_total;
        test;
    end
    methods
        function [E,obj]=init(obj,my_const,electron,E,nm_name)
  
            obj.g_w_period=zeros(1,200);
            obj.test=zeros(1,200);
            obj.g_w_period(1)=electron.g_w;
            obj.g_w_n=1;
            obj.g_w_sum=0;
            %these two comes together to make sure fine and chirp 
            %method gives the same reaults
            obj.z=zeros(size(E.t))';
            obj.N_x=2;
            obj.z(1)=(E.x(2)+0.6*E.x(3))/2;
            obj.x_range=0;
            obj.c=my_const.c;
            obj.lambda_0=my_const.lambda_0;
            obj.dt=E.t(2)-E.t(1);

            
            
           
            obj.v=zeros(size(E.t));
            obj.g_h=electron.g_h;
            obj.gap_h=electron.gap_h;
            obj.n_silica=electron.n_silica;

             
            obj.p_array=zeros(size(E.t));
            gamma=1/sqrt(1-electron.v_0^2/E.c^2);
            obj.W_total=0;
            
            switch nm_name

               case {'tpf_stc_chirp','tpf_standard_chirp','direct'}
                        obj.iteration=@Euler_chirp;
                        k_z_decay=sqrt((2*pi/obj.g_w_period(1))^2-(2*pi/obj.lambda_0).^2);
                         A_evan=exp(-k_z_decay*obj.gap_h/2);     
                        obj.field_inst=A_evan.*real(E.field(E,obj.z(1),E.t(1)));
                        obj.v(1)=electron.v_0-0.5*obj.dt*E.e*obj.field_inst/E.m/gamma^3;
            % p is defined on the same mesh as v on 1+1/2,2+1/2
                        obj.p_array(1)=E.m*obj.v(1)/sqrt(1-obj.v(1)^2/E.c^2);
                        E.E_z_electron(1)=real(obj.field_inst);
            
              
                case 'Euler_matched'
                      %the perfect matched e field;
                      %every position is the negative E field 
                       obj.iteration=@Euler_matched;
                       obj.field_inst=-abs(E.field(E,obj.z(1),E.t(1)));
                       obj.v(1)=electron.v_0-0.5*obj.dt*E.e*obj.field_inst/E.m/gamma^3;
                       % p is defined on the same mesh as v on 1+1/2,2+1/2
                       obj.p_array(1)=E.m*obj.v(1)/sqrt(1-obj.v(1)^2/E.c^2);
                       E.E_z_electron(1)=real(obj.field_inst);
       
        
            end
        end
        

        
        [E,obj]=Euler_chirp(obj,nn,E)

        
    end

end
