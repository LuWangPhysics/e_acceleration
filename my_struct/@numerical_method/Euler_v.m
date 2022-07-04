        function obj=Euler_v(obj,E,nn,field_inst)
                %update p, because it become an explicit method, if v is
                %updated, the right hand side is a function of v too, the
                %method becomes implicit
               obj.p_array(nn)=obj.p_array(nn-1)-obj.dt*E.e*field_inst;
               obj.v(nn)=obj.p_array(nn)/E.m/sqrt(1+obj.p_array(nn)^2/E.c^2/E.m^2);

        end