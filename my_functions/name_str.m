function save_name=name_str(name_head,case_name,electron,E,my_cons,energy)
save_name=['my_output/' case_name name_head num2str(electron.energy_e/electron.e/1e3) 'keV_GDD=' ...
    num2str(my_cons.GDD*1e26) 'e-26_TOD=' num2str(my_cons.TOD*1e40) 'e-40_'   ...
   'E0=' num2str(E.E0*1e-9) 'GVperm_energy' num2str(energy*1e3) 'mJ'];
end