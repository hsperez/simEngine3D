clc;
clear all;

tic

Length = 4;
Width = 0.05;
Volume = 0.05*0.05*Length;
Density = 7800;
Mass = Volume*Density;

gravity = -9.81;

a = Length;
b = Width;
c = Width; 

M = Mass*eye(3);
J_bar = zeros(3,3);

J_bar(1,1) = Mass*(b^2+c^2)/12;
J_bar(2,2) = Mass*(a^2+c^2)/12;
J_bar(3,3) = Mass*(b^2+a^2)/12;


%%%%%  WHEN BODY i is the ground
a_i_bar_x = [1; 0; 0];
a_i_bar_y = [0; 1; 0];
a_i_bar_z = [0; 0; 1];

r_i_o = [0; 0; 0];

p_i_o = [1; 0; 0; 0];
p_i_dot = [0; 0; 0; 0];

s_i_bar = [0; 0; 0];

%%%%%%  INITIAL GUESSES for r_j and p_j
xo = 0;
yo = -1;
zo = -1;

r_j_o = [0; sqrt(2); -sqrt(2)];

p_j_o = [0.653281479071633;
         0.270598051467572;
         0.653281485804746;
         0.270598048678625];

a_j_bar_x = [1; 0; 0];
a_j_bar_y = [0; 1; 0];
a_j_bar_z = [0; 0; 1];

L = 2; 

s_j_bar = [-L;0 ; 0];

c_x = [1; 0; 0];
c_y = [0; 1; 0];
c_z = [0; 0; 1];

r_i = r_i_o;

p_i = p_i_o;


Total_t = 10;
steps = 1000;

h = Total_t/steps;

q_j = zeros(7,steps);
q_j_dot = zeros(7,steps);
q_j_ddot = zeros(7,steps);

Lambdas = zeros(7,steps);

q_j(:,1) = [r_j_o; p_j_o];

F_j = [0; 0; Mass*gravity];

big = zeros(14,14);

%%
for k =  1:steps
    
    i = 1;
    
    error = 1;
    
    t = h*(k-1);
    
    [f, df, ddf] = fun(t);
    
    if k == 1
        
       CD_1 =  GconCD(c_x, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
       CD_2 =  GconCD(c_y, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
       CD_3 =  GconCD(c_z, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
       DP1_1 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_z,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
       DP1_2 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
       DP1_3 = GconDP1(f,df,ddf,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_x,zeros(4,1), q_j_dot(4:7,k));
       e_con = eGcon(q_j(4:7,k),q_j_dot(4:7,k));
       
       Jacobian = [CD_1.Phi_r_j,CD_1.Phi_p_j;
                   CD_2.Phi_r_j,CD_2.Phi_p_j;
                   CD_3.Phi_r_j,CD_3.Phi_p_j;
                   DP1_1.Phi_r_j,DP1_1.Phi_p_j;
                   DP1_2.Phi_r_j,DP1_2.Phi_p_j;
                   DP1_3.Phi_r_j,DP1_3.Phi_p_j;
                   e_con.Phi_r_j, e_con.Phi_p_j];
       
       Phi = [CD_1.Phi;
              CD_2.Phi;
              CD_3.Phi;
              DP1_1.Phi;
              DP1_2.Phi;
              DP1_3.Phi];
          
       Phi_r = Jacobian(1:6,1:3);
       Phi_p = Jacobian(1:6,4:7);
       Phi_econ = Jacobian(7,:);
        
       G = get_G(q_j(4:7,k));
       G_dot = get_G(q_j_dot(4:7,k));
       
       J_P = 4*(G.')*J_bar*G;
       P = q_j(4:7,k).';
       
    
       tau_hat_j = 8*G_dot.'*J_bar*G_dot*q_j(4:7,k);
       gamma_blue = [e_con.gamma;
                     CD_1.gamma;
                     CD_2.gamma;
                     CD_3.gamma;
                     DP1_1.gamma;
                     DP1_2.gamma;
                     DP1_3.gamma];
                       
       
       right_side_blue = [F_j;
                          tau_hat_j;
                          gamma_blue];
       
       
       big(1:3,:) = [M ,zeros(3,4) ,zeros(3,1) ,Phi_r.'] ;
       big(4:7,:) = [zeros(4,3) ,(J_P) ,P.' ,Phi_p.'];
       big(8,:) = [zeros(1,3), P, zeros(1,7)];
       big(9:14,:) = [Phi_r, Phi_p, zeros(6,7)];
        
       red_side = inv(big)*right_side_blue;
       
       q_j_ddot(:,k) = red_side(1:7,1);
       Lambdas(:,k) = red_side(8:14,1);
       
       C_r_n = q_j(1:3,k)+h*q_j_dot(1:3,k);
       C_p_n = q_j(4:7,k)+h*q_j_dot(4:7,k);
       C_r_n_dot = q_j_dot(1:3,k);
       C_p_n_dot = q_j_dot(4:7,k);
       
    
    elseif k ==2 
        
        q_j_ddot(:,k) = q_j_ddot(:,k-1); 
        Lambdas(:,k) = Lambdas(:,k-1);
        
        
        
               %%%%%    BDF     %%%%%%%
            q_j(1:3,k) = C_r_n+(h^2)*q_j_ddot(1:3,k);
            q_j(4:7,k) = C_p_n+(h^2)*q_j_ddot(4:7,k);
            q_j_dot(1:3,k) = C_r_n_dot+h*q_j_ddot(1:3,k);
            q_j_dot(4:7,k) = C_p_n_dot+h*q_j_ddot(4:7,k);
            
                   %%%%%  Residual  %%%%%%%
        while error > 0.00001 
            CD_1 =  GconCD(c_x, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            CD_2 =  GconCD(c_y, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            CD_3 =  GconCD(c_z, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            DP1_1 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_z,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
            DP1_2 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
            DP1_3 = GconDP1(f,df,ddf,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_x,zeros(4,1), q_j_dot(4:7,k));
            e_con = eGcon(q_j(4:7,k),q_j_dot(4:7,k));
            
            Jacobian = [CD_1.Phi_r_j,CD_1.Phi_p_j;
                   CD_2.Phi_r_j,CD_2.Phi_p_j;
                   CD_3.Phi_r_j,CD_3.Phi_p_j;
                   DP1_1.Phi_r_j,DP1_1.Phi_p_j;
                   DP1_2.Phi_r_j,DP1_2.Phi_p_j;
                   DP1_3.Phi_r_j,DP1_3.Phi_p_j;
                   e_con.Phi_r_j, e_con.Phi_p_j];
       
           Phi = [CD_1.Phi;
                  CD_2.Phi;
                  CD_3.Phi;
                  DP1_1.Phi;
                  DP1_2.Phi;
                  DP1_3.Phi];

           Phi_r = Jacobian(1:6,1:3);
           Phi_p = Jacobian(1:6,4:7);
           Phi_econ = e_con.Phi;

           G = get_G(q_j(4:7,k));
           G_dot = get_G(q_j_dot(4:7,k));

           J_P = 4*(G.')*J_bar*G;
           P = q_j(4:7,k).';


           tau_hat_j = 8*G_dot.'*J_bar*G_dot*q_j(4:7,k);
           gamma_blue = [e_con.gamma;
                         CD_1.gamma;
                         CD_2.gamma;
                         CD_3.gamma;
                         DP1_1.gamma;
                         DP1_2.gamma;
                         DP1_3.gamma];


           right_side_blue = [F_j;
                              tau_hat_j;
                              gamma_blue];


           big(1:3,:) = [M ,zeros(3,4) ,zeros(3,1) ,Phi_r.'] ;
           big(4:7,:) = [zeros(4,3) ,(J_P),P.' ,Phi_p.'];
           big(8,:) = [zeros(1,3), P, zeros(1,7)];
           big(9:14,:) = [Phi_r, Phi_p, zeros(6,7)];
           
              
     
           g_residual = [M*q_j_ddot(1:3,k)+(Phi_r.')*Lambdas(2:7,1)-F_j;
                         J_P*q_j_ddot(4:7,k)+(Phi_p.')*Lambdas(2:7,k)+(P.')*Lambdas(1,k)-tau_hat_j;
                         (1/h^2)*Phi_econ;
                         (1/h^2)*Phi];
                         
                         
             test1 = -inv(big);
                     
            correction = -inv(big)*g_residual;
            
            q_j_ddot(:,k) = q_j_ddot(:,k)+correction(1:7,1);
            Lambdas(:,k) = Lambdas(:,k)+correction(8:14,1);
            
            error = norm(correction);
            
            
                     
            i = i+1;

            if i==10  
                error = 0;
            end
            
           end
        
       C_r_n = (4/3)*q_j(1:3,k)-(1/3)*q_j(1:3,k-1)+(8/9)*h*q_j_dot(1:3,k)-(2/9)*h*q_j_dot(1:3,k-1);
       C_p_n = (4/3)*q_j(4:7,k)-(1/3)*q_j(4:7,k-1)+(8/9)*h*q_j_dot(4:7,k)-(2/9)*h*q_j_dot(4:7,k-1);
       C_r_n_dot = (4/3)*q_j_dot(1:3,k)-(1/3)*q_j_dot(1:3,k-1);
       C_p_n_dot = (4/3)*q_j_dot(4:7,k)-(1/3)*q_j_dot(4:7,k-1);
        
    else 
        
        
        q_j_ddot(:,k) = q_j_ddot(:,k-1); 
        Lambdas(:,k) = Lambdas(:,k-1);
        
        
        
               %%%%%    BDF     %%%%%%%
            q_j(1:3,k) = C_r_n+(h^2)*q_j_ddot(1:3,k);
            q_j(4:7,k) = C_p_n+(h^2)*q_j_ddot(4:7,k);
            q_j_dot(1:3,k) = C_r_n_dot+h*q_j_ddot(1:3,k);
            q_j_dot(4:7,k) = C_p_n_dot+h*q_j_ddot(4:7,k);
            
                   %%%%%  Residual  %%%%%%%
        while error > 0.00001 
            CD_1 =  GconCD(c_x, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            CD_2 =  GconCD(c_y, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            CD_3 =  GconCD(c_z, 0, 0, 0, r_i, q_j(1:3,k), p_i, q_j(4:7,k), s_i_bar, s_j_bar,p_i_dot, q_j_dot(4:7,k));
            DP1_1 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_z,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
            DP1_2 = GconDP1(0,0,0,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_z,zeros(4,1), q_j_dot(4:7,k));
            DP1_3 = GconDP1(f,df,ddf,p_i,q_j(4:7,k),a_i_bar_y,a_j_bar_x,zeros(4,1), q_j_dot(4:7,k));
            e_con = eGcon(q_j(4:7,k),q_j_dot(4:7,k));
            
            Jacobian = [CD_1.Phi_r_j,CD_1.Phi_p_j;
                   CD_2.Phi_r_j,CD_2.Phi_p_j;
                   CD_3.Phi_r_j,CD_3.Phi_p_j;
                   DP1_1.Phi_r_j,DP1_1.Phi_p_j;
                   DP1_2.Phi_r_j,DP1_2.Phi_p_j;
                   DP1_3.Phi_r_j,DP1_3.Phi_p_j;
                   e_con.Phi_r_j, e_con.Phi_p_j];
       
           Phi = [CD_1.Phi;
                  CD_2.Phi;
                  CD_3.Phi;
                  DP1_1.Phi;
                  DP1_2.Phi;
                  DP1_3.Phi];

           Phi_r = Jacobian(1:6,1:3);
           Phi_p = Jacobian(1:6,4:7);
           Phi_econ = e_con.Phi;

           G = get_G(q_j(4:7,k));
           G_dot = get_G(q_j_dot(4:7,k));

           J_P = 4*(G.')*J_bar*G;
           P = q_j(4:7,k).';


           tau_hat_j = 8*G_dot.'*J_bar*G_dot*q_j(4:7,k);
           gamma_blue = [e_con.gamma;
                         CD_1.gamma;
                         CD_2.gamma;
                         CD_3.gamma;
                         DP1_1.gamma;
                         DP1_2.gamma;
                         DP1_3.gamma];


           right_side_blue = [F_j;
                              tau_hat_j;
                              gamma_blue];


           big(1:3,:) = [M ,zeros(3,4) ,zeros(3,1) ,Phi_r.'] ;
           big(4:7,:) = [zeros(4,3) ,(J_P),P.' ,Phi_p.'];
           big(8,:) = [zeros(1,3), P, zeros(1,7)];
           big(9:14,:) = [Phi_r, Phi_p, zeros(6,7)];
           
              
     
           g_residual = [M*q_j_ddot(1:3,k)+(Phi_r.')*Lambdas(2:7,1)-F_j;
                         J_P*q_j_ddot(4:7,k)+(Phi_p.')*Lambdas(2:7,k)+(P.')*Lambdas(1,k)-tau_hat_j;
                         (1/h^2)*Phi_econ;
                         (1/h^2)*Phi];
                         
                         
             test1 = -inv(big);
                     
            correction = -inv(big)*g_residual;
            
            q_j_ddot(:,k) = q_j_ddot(:,k)+correction(1:7,1);
            Lambdas(:,k) = Lambdas(:,k)+correction(8:14,1);
            
            error = norm(correction);
            
            
                     
            i = i+1;

            if i==10  
                error = 0;
            end
            
           end
        
       C_r_n = (4/3)*q_j(1:3,k)-(1/3)*q_j(1:3,k-1)+(8/9)*h*q_j_dot(1:3,k)-(2/9)*h*q_j_dot(1:3,k-1);
       C_p_n = (4/3)*q_j(4:7,k)-(1/3)*q_j(4:7,k-1)+(8/9)*h*q_j_dot(4:7,k)-(2/9)*h*q_j_dot(4:7,k-1);
       C_r_n_dot = (4/3)*q_j_dot(1:3,k)-(1/3)*q_j_dot(1:3,k-1);
       C_p_n_dot = (4/3)*q_j_dot(4:7,k)-(1/3)*q_j_dot(4:7,k-1);
       
        
        
    end
    
    
end
    
    toc
    
    
    




%%





function DP1 = GconDP1(f,df,ddf,p_i,p_j,a_i_bar,a_j_bar,p_i_dot, p_j_dot)

    


    [eo_i, e_i] = e_func(p_i);
    [eo_j, e_j] = e_func(p_j);
    e_i_tilde = tilde(e_i);
    e_j_tilde = tilde(e_j);

    A_i = eye(3)*(2*(eo_i^2)-1)+2*(e_i*e_i.'+eo_i*e_i_tilde);
    A_j = eye(3)*(2*(eo_j^2)-1)+2*(e_j*e_j.'+eo_j*e_j_tilde);
    
    a_i = A_i*a_i_bar;
    a_j = A_j*a_j_bar;
    
    
    B_i = B(p_i,a_i_bar);
    B_j = B(p_j,a_j_bar);
   
    
    B_i_dot = B(p_i_dot,a_i_bar);
    B_j_dot = B(p_j_dot,a_j_bar);
    
    a_i_dot = B_i*p_i_dot;
    a_j_dot = B_j*p_j_dot;
    
    B_i_bar = B(p_i,a_i_bar);
    B_j_bar = B(p_j,a_j_bar);
    
    
    DP1.Phi = a_i.'*a_j - f;
    
    DP1.nu = df;
    
    DP1.gamma = -(2*a_i_dot.'*a_j_dot)-(a_j.'*B_i_dot*p_i_dot)-(a_i.'*B_j_dot*p_j_dot)+ddf;
    
    
    DP1.Phi_r_i = [0 0 0];
    DP1.Phi_r_j = [0 0 0];
    
    DP1.Phi_r = [DP1.Phi_r_i, DP1.Phi_r_j];
    
    DP1.Phi_p_i = a_j.'*B_i_bar;
    DP1.Phi_p_j = a_i.'*B_j_bar;
    
    DP1.Phi_p = [DP1.Phi_p_i,DP1.Phi_p_j];


end 

function CD = GconCD(c, f, df, ddf, r_i, r_j, p_i, p_j,s_i_P_bar, s_j_Q_bar, p_i_dot, p_j_dot)
                    
    [eo_i,e_i] = e_func(p_i);
    [eo_j,e_j] = e_func(p_j);
    e_i_tilde = tilde(e_i);
    e_j_tilde = tilde(e_j);
    
    A_i = eye(3)*(2*(eo_i^2)-1)+2*(e_i*e_i.'+eo_i*e_i_tilde);
    A_j = eye(3)*(2*(eo_j^2)-1)+2*(e_j*e_j.'+eo_j*e_j_tilde);
    
    s_i_P = A_i*s_i_P_bar;
    s_j_Q = A_j*s_j_Q_bar;
    
   
    d_ij = getdij(r_i,r_j,s_i_P,s_j_Q);
   
    
     B_j_dot = B(p_j_dot,s_j_Q_bar);
     B_i_dot = B(p_i_dot,s_i_P_bar);
    
   

    CD.Phi = c.'*d_ij -f;
    CD.nu = df;
    CD.gamma = -c.'*B_j_dot*p_j_dot+c.'*B_i_dot*p_i_dot+ddf;
    CD.Phi_r_i = -c.';
    CD.Phi_r_j = c.';
    CD.Phi_r = [CD.Phi_r_i, CD.Phi_r_j];
    CD.Phi_p_i = -c.'*B(p_i,s_i_P_bar);
    CD.Phi_p_j = c.'*B(p_j,s_j_Q_bar);
    CD.Phi_p = [CD.Phi_p_i,CD.Phi_p_j];
    


end

function e_con = eGcon(p_i, p_i_dot)

    e_con.Phi = p_i.'*p_i-1;

    e_con.Phi_r_j = zeros(1,3);

    e_con.Phi_p_j = 2*(p_i.');
    
    e_con.Phi_r_p = [e_con.Phi_r_j, e_con.Phi_p_j];
    
    e_con.nu = 0;
    
    e_con.gamma = -2*(p_i_dot.')*p_i_dot;

end


function dij = getdij(r_i,r_j,s_i_P,s_j_P)

    dij = r_j+s_j_P- r_i -s_i_P;
    
end

function Bop = B(p,a)

    eo = p(1);
    e = [p(2); p(3); p(4)];
    atilde = tilde(a);
    etilde = tilde(e);

    Bop = 2*[(eo*eye(3)+etilde)*a, e*(a.')-(eo*eye(3)+etilde)*atilde];
    
end


function [eo, e] = e_func(f)

    eo = f(1);
    e = [f(2); f(3); f(4)];
    
end

function Matrix_t = tilde(f)

Matrix_t = [0 -f(3) f(2); f(3) 0 -f(1); -f(2) f(1) 0];

end

function [f, df, ddf] = fun(t)

    f = sin((pi/4)*cos(2*t));
    df = -pi*sin(2*t)*cos(pi*cos(2*t)/4)/2;
    ddf = -pi*cos(2*t)*cos(pi*cos(2*t)/4) - pi^2*sin(2*t)^2*sin(pi*cos(2*t)/4)/4;

end

function G = get_G(p)

    G = [
        -p(2) p(1) p(4) -p(3);
        -p(3) -p(4) p(1) p(2);
        -p(4) p(3) -p(2)  p(1)];

end






