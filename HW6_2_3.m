clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                     INPUTS
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%

DP1results = GconDP1(f,df,ddf,p_i,p_j,a_i_bar,a_j_bar,p_i_dot, p_j_dot)

CDresults = GconCD(c, f, df, ddf, r_i, r_j, p_i, p_j, s_i_P_bar, s_j_Q_bar,p_i_dot, p_j_dot)



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
    
    
    Phi_r_i = [0 0 0];
    Phi_r_j = [0 0 0];
    
    DP1.Phi_r = [Phi_r_i, Phi_r_j];
    
    Phi_p_i = a_j.'*B_i_bar;
    Phi_p_j = a_i.'*B_j_bar;
    
    DP1.Phi_p = [Phi_p_i,Phi_p_j];


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
    Phi_r_i = -c.';
    Phi_r_j = c.';
    CD.Phi_r = [Phi_r_i, Phi_r_j];
    Phi_p_i = -c.'*B(p_i,s_i_P_bar);
    Phi_p_j = c.'*B(p_j,s_j_Q_bar);
    CD.Phi_p = [Phi_p_i,Phi_p_j];
    


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


