function [fwh_m]=funfwhm(t,Et)
position=find(abs(Et)==max(abs(Et))); %find position of the peak
peak=Et(position);
p_position=t(position);
position_2=find((abs(Et)>max(abs(Et)/2)));%find range of the position that is within fwhm
fwh_m=abs(t(position_2(1))-t(position_2(end)));
end