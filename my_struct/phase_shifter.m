classdef phase_shifter
    properties
        shifter_n;
        phase_shifter_location;
        shifter_flag;
        c,g_w;
        ratio;
    end
    methods
        
        function obj=init(obj,N,my_cons,electron,E0)
            obj.c=my_cons.c;
            obj.g_w=electron.g_w;
            obj.phase_shifter_location=zeros(1,N);
            obj.phase_shifter_location(1)=1;
            obj.shifter_n=1;
            obj.shifter_flag=1;
            delta_v=abs(E0*electron.e/electron.m)*obj.g_w/electron.v_0/10;
            if electron.energy_e<1e3*electron.e
                obj.ratio=0.98;
            else
            obj.ratio=1-(delta_v)/electron.v_0;
            end
        end
        
        function obj=flag_calculate(obj,NM,nn)
               if NM.v(nn-1)<(max(NM.v)*abs(obj.ratio))
      
                       loc_id=max(obj.phase_shifter_location);
                       %the distance electron travels, including backwards
                       z_path=sum(abs(diff(NM.z(loc_id:(nn-1)))));
                            if (z_path>2*obj.g_w)||(NM.v(nn-1)<0)          
                               obj.shifter_n=obj.shifter_n+1;
                               obj.phase_shifter_location(obj.shifter_n)=nn-1;
                               obj.shifter_flag=obj.shifter_flag;
                            end
               end
        end
        
        
        
    end
end