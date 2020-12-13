
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                    FLUID FORCE                          %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef    Flf 
    methods     ( Static = true )
        
        function    F = get_FlF(q_i,q_i_dot,Nsph,R,nu,sr,SA)
           
            F = zeros(Nsph*3,1);
            
            for i = 1:Nsph
                
                x = 3*(i-1)+1;
                y = 3*(i-1)+2;
                z = 3*(i-1)+3;
                
                
                F(x,1) = get_xF(q_i(y),sr,q_i_dot(x),R,nu);
                F(y,1) = get_yF(q_i(y),sr,q_i_dot(y),R,nu,SA);
                F(z,1) = 0;
                

                
            end
        end
        

    end
end

function xF = get_xF(y,sr,vx,R,nu)

    xF = 6*pi*nu*R*(y*sr-vx);

end

function yF = get_yF(y,sr,vy,R,nu,SA)

    yF = nu*SA*(-abs((y+sqrt(3)*R/2)*sr/(2*R))+abs((y-sqrt(3)*R/2)*sr)/(2*R))-6*pi*nu*R*(vy);

end


