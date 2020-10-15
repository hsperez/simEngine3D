clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                     INPUTS
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sym t;
% 
% theta = (pi/4)cos(2*t);

%%%%%  WHEN BODY i is the ground
a_i_bar_x = [1; 0; 0];
a_i_bar_y = [0; 1; 0];
a_i_bar_z = [0; 0; 1];

r_i_o = [0; 0; 0];

p_i_o = [1; 0; 0; 0];

s_i_bar = [0; 0; 0];

%%%%%%  INITIAL GUESSES for r_j and p_j
xo = 0;
yo = -1;
zo = -1;

r_j_o = [xo; yo; zo];

p_j_o = [0.707; -.2; .7; -.2];

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

%%%% Setting newton Raphson
n = 100;
i = 2;

r_j = zeros(3,n);
p_j = zeros(4,n);
r_p_j = zeros(7,n);

r_j(:,1) = r_j_o;
p_j(:,1) = p_j_o;

r_p_j(:,1) = [r_j(:,1).',p_j(:,1).'].';



    t = 0; 

    [f, df, ddf] = fun(t);


    error = 1;
    error_r = zeros(7,1);

while error > 0.000001

   CD_1 =  GconCD(c_x, 0, 0, 0, r_i, r_j(:,i-1), p_i, p_j(:,i-1), s_i_bar, s_j_bar);
   CD_2 =  GconCD(c_y, 0, 0, 0, r_i, r_j(:,i-1), p_i, p_j(:,i-1), s_i_bar, s_j_bar);
   CD_3 =  GconCD(c_z, 0, 0, 0, r_i, r_j(:,i-1), p_i, p_j(:,i-1), s_i_bar, s_j_bar);
   DP1_1 = GconDP1(0,0,0,p_i,p_j(:,i-1),a_i_bar_z,a_j_bar_z);
   DP1_2 = GconDP1(0,0,0,p_i,p_j(:,i-1),a_i_bar_y,a_j_bar_z);
   DP1_3 = GconDP1(sin(f),0,0,p_i,p_j(:,i-1),a_i_bar_y,a_j_bar_x);
   e_con = eGcon(p_j(:,i-1));
   
   
   Multiplier = [r_p_j(1,i-1) 0 0 0 0 0 0;
                 0 r_p_j(2,i-1) 0 0 0 0 0;
                 0 0 r_p_j(3,i-1) 0 0 0 0;
                 0 0 0 r_p_j(4,i-1) 0 0 0;
                 0 0 0 0 r_p_j(5,i-1) 0 0;
                 0 0 0 0 0 r_p_j(6,i-1) 0;
                 0 0 0 0 0 0 r_p_j(7,i-1)];
   
   
   Jacobian = [CD_1.Phi_r_j,CD_1.Phi_p_j;
               CD_2.Phi_r_j,CD_2.Phi_p_j;
               CD_3.Phi_r_j,CD_3.Phi_p_j;
               DP1_1.Phi_r_j,DP1_1.Phi_p_j;
               DP1_2.Phi_r_j,DP1_2.Phi_p_j;
               DP1_3.Phi_r_j,DP1_3.Phi_p_j;
               e_con.Phi_r_j, e_con.Phi_p_j];
           
           
    Ultimate_Phi = [CD_1.Phi;
                    CD_2.Phi;
                    CD_3.Phi;
                    DP1_1.Phi;
                    DP1_2.Phi;
                    DP1_3.Phi;
                    e_con.Phi];
    
    
    
    r_p_j(:,i) = r_p_j(:,i-1)-(inv(Jacobian))*Ultimate_Phi;
    
    r_j(:,i) = r_p_j(1:3,i);
    p_j(:,i) = r_p_j(4:7,i);
    
    
    
   for m = 1:7
   
       error_r(m) = (r_p_j(m,i)-r_p_j(m,i-1))/(r_p_j(m,i));
       
   end
   
   error = max(error_r);
    
   i = i +1;
   
    if i== 100
        error = 0;
    end 
end

r_j = r_j(:,1:(i-1));
p_j = p_j(:,1:(i-1));
r_p_j = r_p_j(:,1:(i-1));

%%

%DP1 = GconDP1(f,df,ddf,p_i,p_j,a_i_bar,a_j_bar,p_i_dot, p_j_dot)

%CD = GconCD(c, f, df, ddf, r_i, r_j, p_i, p_j, s_i_P_bar, s_j_Q_bar,p_i_dot, p_j_dot)



%%

function DP1 = GconDP1(f,df,ddf,p_i,p_j,a_i_bar,a_j_bar) %,p_i_dot, p_j_dot)

    


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
   
    
%     B_i_dot = B(p_i_dot,a_i_bar);
%     B_j_dot = B(p_j_dot,a_j_bar);
    
%     a_i_dot = B_i*p_i_dot;
%     a_j_dot = B_j*p_j_dot;
    
    B_i_bar = B(p_i,a_i_bar);
    B_j_bar = B(p_j,a_j_bar);
    
    
    DP1.Phi = a_i.'*a_j - f;
    
    DP1.nu = df;
    
    %DP1.gamma = -(2*a_i_dot.'*a_j_dot)-(a_j.'*B_i_dot*p_i_dot)-(a_i.'*B_j_dot*p_j_dot)+ddf;
    
    
    DP1.Phi_r_i = [0 0 0];
    DP1.Phi_r_j = [0 0 0];
    
    DP1.Phi_r = [DP1.Phi_r_i, DP1.Phi_r_j];
    
    DP1.Phi_p_i = a_j.'*B_i_bar;
    DP1.Phi_p_j = a_i.'*B_j_bar;
    
    DP1.Phi_p = [DP1.Phi_p_i,DP1.Phi_p_j];


end 

function CD = GconCD(c, f, df, ddf, r_i, r_j, p_i, p_j,s_i_P_bar, s_j_Q_bar) %, p_i_dot, p_j_dot)
                    
    [eo_i,e_i] = e_func(p_i);
    [eo_j,e_j] = e_func(p_j);
    e_i_tilde = tilde(e_i);
    e_j_tilde = tilde(e_j);
    
    A_i = eye(3)*(2*(eo_i^2)-1)+2*(e_i*e_i.'+eo_i*e_i_tilde);
    A_j = eye(3)*(2*(eo_j^2)-1)+2*(e_j*e_j.'+eo_j*e_j_tilde);
    
    s_i_P = A_i*s_i_P_bar;
    s_j_Q = A_j*s_j_Q_bar;
    
   
    d_ij = getdij(r_i,r_j,s_i_P,s_j_Q);
   
    
%     B_j_dot = B(p_j_dot,s_j_Q_bar);
%     B_i_dot = B(p_i_dot,s_i_P_bar);
    
   

    CD.Phi = c.'*d_ij -f;
    CD.nu = df;
    %CD.gamma = -c.'*B_j_dot*p_j_dot+c.'*B_i_dot*p_i_dot+ddf;
    CD.Phi_r_i = -c.';
    CD.Phi_r_j = c.';
    CD.Phi_r = [CD.Phi_r_i, CD.Phi_r_j];
    CD.Phi_p_i = -c.'*B(p_i,s_i_P_bar);
    CD.Phi_p_j = c.'*B(p_j,s_j_Q_bar);
    CD.Phi_p = [CD.Phi_p_i,CD.Phi_p_j];
    


end

function e_con = eGcon(p_i)

    e_con.Phi = p_i.'*p_i-1;

    e_con.Phi_r_j = zeros(1,3);

    e_con.Phi_p_j = 2*(p_i.');
    
    e_con.Phi_r_p = [e_con.Phi_r_j, e_con.Phi_p_j];

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

    f = (pi/4)*cos(2*t);
    df = (pi/2)*sin(2*t);
    ddf = -pi*cos(2*t);

end


