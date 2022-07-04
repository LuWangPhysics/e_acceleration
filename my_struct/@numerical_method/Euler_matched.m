function [E,obj]=Euler_matched(obj,nn,E)

            
          
                 obj=obj.Euler_z(nn);
                 %calculate force using B and E at mesh 1,2,3,4   
                 obj.field_inst=ps.shifter_flag.*E.field(E,obj.z(nn),E.t(nn));
                 obj.field_inst=-abs(obj.field_inst);
                 %obj.field_inst=E.field(E,obj.z(nn),E.t(nn));
                 obj=obj.Euler_v(E,nn,real(obj.field_inst));
                 E.E_z_electron(nn-1)=obj.field_inst;
            
    
                 
end