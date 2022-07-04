function E_delta_eV =delta_kinetic(NM,my_cons,electron)
            gamma=@(v) 1./sqrt(1-v.^2./my_cons.c^2);
            E_delta_eV=electron.m*my_cons.c^2*(gamma(NM.v)-gamma(NM.v(1)))/electron.e;
end