function [x_range,t_period ,obj]=g_w_period_construct(obj,E,nn,N_margin)
 %if it is the stc, lambda is changing with respect to x
 %get the g_w for the varing lambda
%       if ~isempty(E.d_lambda)
%          dx=E.x(2)-E.x(1);
%         lambda_x= E.d_lambda(fix(obj.z(nn-1)/dx)+1);
%         obj.g_w_period(obj.g_w_n)=lambda_x*(obj.v(nn-1)/obj.c);  
%         obj.g_h=lambda_x/2/(obj.n_silica-1); 
%           
%       else
%       % obj.g_h has the initialization value 
%       obj.g_w_period(obj.g_w_n)=obj.lambda_0*(obj.v(nn-1)/obj.c);                      %spatial grating period
%       end
      
       obj.g_w_period(obj.g_w_n)=obj.lambda_0*(obj.v(nn-1)/obj.c);           

       
       
%if the kinetic energy is too small, set the period to the minimum value     
      if obj.g_w_period(obj.g_w_n)<0.1e-6
        obj.g_w_period(obj.g_w_n)=1e-7;
      end
      
 %--------------------
 %get the delay array
 %--------------------
      
    
      t_silica=obj.g_h*(obj.n_silica-1)/obj.c;
      t_center=(t_silica)/2;
      t_amp=(t_silica)/2;
      dx=E.x(2)-E.x(1);
      shape_sharp=3;
      x=0:dx:(obj.g_w_period(obj.g_w_n)+N_margin*dx);
      x_range=length(x);

      %air to silica time shift
      t_period=t_amp.*tanh(shape_sharp*sin(pi+2*pi.*x./(obj.g_w_period(obj.g_w_n))))+t_center;   
end