        function obj=Euler_z(obj,nn)
                %z,t is defined at mesh 1,2 3...
                %v is defined as 1+1/2, 2+1/2...
                %use z and t to have E, B at 1,2,3,4
                obj.z(nn)=obj.z(nn-1)+obj.v(nn-1)*obj.dt;
 
        end