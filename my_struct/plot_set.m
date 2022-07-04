classdef plot_set
    properties
        name;
        save_flag;
    end
    methods
        function obj=init(obj,name_str,save_flag)
        obj.name=name_str;
        obj.save_flag=save_flag;
        end


        function vz(obj,E,NM,c)
             N_z=find(NM.z==max(NM.z));

            
            figure
            plot(NM.z(1:N_z)./1e-3,real(E.E_z_electron(1:N_z)))
            xlabel('x (mm)')
            ylabel('Field strength (V/m)')
            title('the field electron sees')
                 if obj.save_flag==1
            savefig(gcf,[obj.name 'E_z.fig']);
                 end
                 
                 
                 
        end
        
        function delta_kenetic(obj,NM,my_cons,electron)
            N_z=find(NM.z==max(NM.z));
            gamma=@(v) 1./sqrt(1-v.^2./my_cons.c^2);
            E_delta_eV=electron.m*my_cons.c^2*(gamma(NM.v)-gamma(NM.v(1)))/electron.e;
            figure
             %hold on
            plot(NM.z(1:N_z).*1e3,E_delta_eV(1:N_z)./1e3)
            xlabel('x (mm)')
            ylabel('energy gain (keV)')
            if obj.save_flag==1
               
            savefig(gcf,[obj.name '_d_E.fig']);
            end
        end
        

        
    end
    
end