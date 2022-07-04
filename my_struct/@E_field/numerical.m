function output=numerical(obj,z,t_j)
dt=obj.t(2)-obj.t(1);
%select the postion of x that the z falls in between
temp=obj.x-z;
temp2=zeros(size(temp));
temp2(1:end-1)=temp(2:end);
temp2(end)=temp2(end-1);
t_p=fix((t_j-obj.t(1))./dt)+1;
x_p=find((temp.*temp2)<0);
if isempty(x_p)
    x_p=1;
end
E_before=obj.E_xt(t_p,x_p);
E_after=obj.E_xt(t_p,x_p+1);
dx=obj.x(2)-obj.x(1);
output=(abs(obj.x(x_p)-z)*E_before+abs(obj.x(x_p+1)-z)*E_after)./dx;
end