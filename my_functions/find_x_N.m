function N_x=find_x_N(x_arr,z_electron)
     
       temp=x_arr-z_electron;
      temp2=zeros(size(temp));
      temp2(1:end-1)=temp(2:end);
      temp2(end)=temp2(end-1);

      N_x=find((temp.*temp2)<0)+1;
end