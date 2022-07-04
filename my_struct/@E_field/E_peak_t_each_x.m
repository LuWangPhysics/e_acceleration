function obj=E_peak_t_each_x(obj,z,nn)

%select the postion of x that the z falls in between
temp=obj.x-z;
temp2=zeros(size(temp));
temp2(1:end-1)=temp(2:end);
temp2(end)=temp2(end-1);
N_x_start=0;%NM.N_x-NM.x_range-1;
x_p=find((temp.*temp2)<0);
E_arr=abs(obj.E_xt(:,x_p-N_x_start));
t_p=find(E_arr==max(E_arr));
obj.E_peak_t(nn)=obj.t(t_p(1));
end