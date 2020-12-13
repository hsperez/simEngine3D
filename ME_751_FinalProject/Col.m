
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                         %
%                    COLLISION                            %
%                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef    Col 
    methods     ( Static = true )
        
        function    [F D] = get_LF(q_i,q_i_dot,Nsph,R,nu, rep)
            
            tao = 1/(0.01*R);
            
            A1 = 3*pi*sqrt(2)/4;
            A2 = 231*pi*sqrt(2)/80;
            
            sc = 0.225*R;
            
            F = zeros(Nsph*3,1);
            D_o = zeros(Nsph,1);
            d_o = zeros(Nsph,1);
            
            for i = 1:Nsph
                
                x = 3*(i-1)+1;
                y = 3*(i-1)+2;
                z = 3*(i-1)+3;

                for j = 1:Nsph
                    
                    xj = 3*(j-1)+1;
                    yj = 3*(j-1)+2;
                    zj = 3*(j-1)+3;
                    
                    if i ==j 
                        
                        F(x,1) = F(x,1)+0;
                        F(y,1) = F(y,1)+0;
                        F(z,1) = F(z,1)+0; 
                        
                        d_o(j,1) = 10000;
                        
                    else 
                        
                        d = distance(q_i(x:z),q_i(xj:zj));
                        
                        d_o(j,1) = d;
                        
                        vector = uni_vec(q_i(x:z),q_i(xj:zj),d);
                        vdif = vel_dif(q_i_dot(x:z),q_i_dot(xj:zj),vector);
                        
                        
                        FLub = -(nu/2)*(A1*((2*R/d)^(3/2))+A2*sqrt(2*R/d));
                        
                        if d < sc
                        
                            FLub_c = -(nu/2)*(A1*((2*R/sc)^(3/2))+A2*sqrt(2*R/sc));
                            
                        else 
                            
                            FLub_c = 0;
                        
                        end
                        
                        LForce = (FLub-FLub_c)*vector/1000;
                        
                        if rep == 1
                            
                            RForce = (((FLub-FLub_c)/1000)*tao*exp(-tao*d)/(1-exp(-tao*d)))*vector ;
                            
                        else
                            
                            RForce = 0*vector;
                            
                        end
                        
                        F(x:z,1) = F(x:z,1) + LForce + RForce;
                        
                    end
                end
                
                D_o(i,1) = min(d_o(:,1));
                
            end
            
            D = min(D_o(:,1));
            
        end
        

    end
end

function d = distance(sph1,sph2)

    d = sqrt((sph1(1)-sph2(1))^2+(sph1(2)-sph2(2))^2+(sph1(3)-sph2(3))^2);

end

function vdif = vel_dif(vsph1, vsph2,vector)

    vdif = sqrt((vsph1(1)-vsph2(1))^2+(vsph1(2)-vsph2(2))^2+(vsph1(3)-vsph2(3))^2);
    
    %vdif = (vsph1-vsph2).*(vector.');

end


function vec = uni_vec(sph1,sph2,d)

    vec = [(sph2(1)-sph1(1))/d;
           (sph2(2)-sph1(2))/d;
           (sph2(3)-sph1(3))/d];

end



