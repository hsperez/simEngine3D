clc
clear all
close all

filename = 'Spheres.txt';

% Units are in mm

T = readtable(filename);

Nsph = size(T,1);

Radius = .01; %mm
R = Radius;
Volume = 4*pi*Radius^3/3;
Density = 1;
Mass = Volume*Density;

SA = 4*pi*Radius^2;


%gravity negligeble
sr = 10; %Shear Rate
nu = .112; % kg/(mm*s)

M = Mass*eye(3*Nsph);
J_bar = (Mass*2*Radius^2/5)*eye(3);


a_i_bar_x = [1; 0; 0];
a_i_bar_y = [0; 1; 0];
a_i_bar_z = [0; 0; 1];

r_i = [-.1; .1; 0];
r_j = [.1; -.1; 0];

p_i = [1; 0; 0; 0];
p_j = [1; 0; 0; 0];

s_i_bar = [Radius; 0; 0];


Total_t = 10;
steps = 1000;

h = Total_t/steps;

q_i = zeros(7*Nsph,steps);
q_i_dot = zeros(7*Nsph,steps);
q_i_ddot = zeros(7*Nsph,steps);

F_t = zeros(3,steps);

q_i(:,1) = [r_i; r_j; p_i; p_j];

for k = 1:steps
    
   t = h*(k-1);
   
   i = 1;
   error = 1;
   
   if k == 1
       
        G_i = get_G(q_i(7:10,k));
        G_j = get_G(q_i(11:14,k));
        G_dot_i = get_G(q_i(7:10,k));
        G_dot_j = get_G(q_i(11:14,k));
        
        J_P_i = 4*(G_i.')*J_bar*G_i;
        J_P_j = 4*(G_j.')*J_bar*G_j;
       
        
        tau_hat_i = 8*G_dot_i.'*J_bar*G_dot_i*q_i(7:10,k);
        tau_hat_j = 8*G_dot_j.'*J_bar*G_dot_j*q_i(11:14,k);
        
        tau_hat = [tau_hat_i;
                   tau_hat_j];
               
              
        
        xF_i = get_xF(q_i(2,k),sr,q_i_dot(1,k),R,nu,SA);
        yF_i = get_yF(q_i(2,k),sr,q_i_dot(2,k),R,nu,SA);
        
        xF_j = get_xF(q_i(5,k),sr,q_i_dot(4,k),R,nu,SA);
        yF_j = get_yF(q_i(5,k),sr,q_i_dot(5,k),R,nu,SA);
         
        F_i = [xF_i; yF_i; 0];
        F_j = [xF_j; yF_j; 0];
        
        F_t(:,k) = F_i;
        F = [F_i;F_j];
        
        right_side_blue = [F;
                          tau_hat];
                      
        big = [M, zeros(6,8);
               zeros(4,6), J_P_i,zeros(4,4);
               zeros(4,6), zeros(4,4),J_P_j];
           
        red_side = right_side_blue\big;
        
        q_i_ddot(:,k) = red_side;
        
        C_r_n = q_i(1:6,k)+h*q_i_dot(1:6,k);
        C_p_n = q_i(7:14,k)+h*q_i_dot(7:14,k);
        C_r_n_dot = q_i_dot(1:6,k);
        C_p_n_dot = q_i_dot(7:14,k);
        
       
   else
       
        q_i_ddot(:,k) = q_i_ddot(:,k-1);
        
        %%%%%%%%%%%%    BDF    %%%%%%%%%%%%%%

        q_i(1:6,k)= C_r_n+(h^2)*q_i_ddot(1:6,k);
        q_i(7:14,k) = C_p_n+(h^2)*q_i_ddot(7:14,k);
        q_i_dot(1:6,k) = C_r_n_dot+h*q_i_ddot(1:6,k);
        q_i_dot(7:14,k) = C_p_n_dot+h*q_i_ddot(7:14,k);
        
        %%%%%%%%%%%%    RESIDUAL    %%%%%%%%%%%%%%
        
        while error > 0.001 
        
        G_i = get_G(q_i(7:10,k));
        G_j = get_G(q_i(11:14,k));
        G_dot_i = get_G(q_i(7:10,k));
        G_dot_j = get_G(q_i(11:14,k));
        
        J_P_i = 4*(G_i.')*J_bar*G_i;
        J_P_j = 4*(G_j.')*J_bar*G_j;
        
        J_P = [J_P_i, zeros(4,4);
               zeros(4,4), J_P_j];
       
        tau_hat_i = 8*G_dot_i.'*J_bar*G_dot_i*q_i(7:10,k);
        tau_hat_j = 8*G_dot_j.'*J_bar*G_dot_j*q_i(11:14,k);
        
        tau_hat = [tau_hat_i;
                   tau_hat_j];
        
        xF_i = get_xF(q_i(2,k),sr,q_i_dot(1,k),R,nu,SA);
        yF_i = get_yF(q_i(2,k),sr,q_i_dot(2,k),R,nu,SA);
        
        xF_j = get_xF(q_i(5,k),sr,q_i_dot(4,k),R,nu,SA);
        yF_j = get_yF(q_i(5,k),sr,q_i_dot(5,k),R,nu,SA);
         
        F_i = [xF_i; yF_i; 0];
        F_j = [xF_j; yF_j; 0];
        
        F = [F_i;F_j];
        
        F_t(:,k) = F_i;
        
        right_side_blue = [F;
                          tau_hat];
                      
        big = [M, zeros(3*Nsph,4*Nsph);
               zeros(4,6), J_P_i,zeros(4,4);
               zeros(4,6), zeros(4,4),J_P_j];
           
        g_residual = [M*q_i_ddot(1:6,k)-F;
                      J_P*q_i_ddot(7:14,k)-tau_hat];
           
        red_side = right_side_blue\big;
        
        q_i_ddot(:,k) = red_side;
        
        correction = -g_residual\big;
        
        error = norm(correction);
        
        i = i+1;
        
            if i==100
                error = 0;

            end 
        end

        
        C_r_n = q_i(1:6,k)+h*q_i_dot(1:6,k);
        C_p_n = q_i(7:14,k)+h*q_i_dot(7:14,k);
        C_r_n_dot = q_i_dot(1:6,k);
        C_p_n_dot = q_i_dot(7:14,k);
        
        

       
   end
       
   
   
    
end


function xF = get_xF(y,sr,vx,R,nu,SA)

    xF = 6*pi*nu*R*(y*sr-vx);

end

function yF = get_yF(y,sr,vx,R,nu,SA)

    yF = nu*SA*(-abs((y+sqrt(3)*R/2)*sr/(2*R))+abs((y-sqrt(3)*R/2)*sr)/(2*R));

end


function G = get_G(p)

    G = [
        -p(2) p(1) p(4) -p(3);
        -p(3) -p(4) p(1) p(2);
        -p(4) p(3) -p(2)  p(1)];

end

