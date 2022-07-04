function [obj,x_start,x_range]=e_x_select(obj,x_arr,e_xt_out)
          dx=x_arr(2)-x_arr(1);
          %shift t postion
          N_center=fix(length(x_arr)/2);
  

          x_shift=obj.z_center; 
          x_start=N_center-fix(x_shift/dx)-1;
          x_temp=x_arr+x_shift;
          obj.x=x_temp(x_start:end);
          x_range=1:length(obj.x);
          obj.t_start_px=find(abs(e_xt_out(:,x_start))==max(abs(e_xt_out(:,x_start))));

          
end