function TF_with_feedback_LQR
    % Define parameters
    c = 0.5;   % Damping coefficient
    m = 0.8;   % Mass
    g = 9.81;  % Gravitational acceleration
    dc = 1;  % Distance from pivot to center of mass
    d1 = 2;    % Distance from pivot to thrust point
    m_motor = 0.2; % Mass of motor
    J = 1/3*(m*d1^2) + m_motor*d1^2; % Moment of inertia
    
    % Define the linearized system matrices A and B
    A = [0, 1; -m*g*dc/J, -c/J];
    B = [0; d1/J];
    
    % Output matrix C
    C = [1, 0];
    
    % D matrix 
    D = 0;
    
   
    % Transfer function computation
    sys = ss(A, B, C, D);
    disp(A);
    disp(B);
    transferFunction = tf(sys);

    co = ctrb(sys);
    controllability = rank(co); 
    
    disp("controllability matrix : ");
    disp(co); 
    
    [V,Diag] = eig(co);
    disp("eigenvector of co");
    disp(V);

    disp("eigenvalues of co");
    disp(Diag);

    disp("Controllabilite : ");
    disp(controllability);

    Wc = gram(sys,'c');
    disp("Gramian : ");
    disp(Wc);
    disp(eig(Wc));

    [U,S,V] = svd(co,'econ');
    disp("SVD U : ");
    disp(U);

    disp("SVD S : ");
    disp(S);
    
    disp("SVD V : ");
    disp(V);
    
    % Design LQR controller
    Q = [10 0; 0 1];

    disp(Q);

    R = 1;     
    K = lqr(A, B, Q, R);

    % Display the LQR controller gains
    disp('LQR Controller Gains:');
    disp(K);

    disp(Q);

    Ac = [(A-B*K)];
    Bc = [B];
    Cc= [C];
    Dc = [D];
    
    disp(eig(A));
    disp(eig(Ac));
    figure
    bode(sys);
    sys_cl = ss(Ac,Bc,Cc,Dc);

    [Gm_lqr, Pm_lqr, Wcg_lqr, Wcp_lqr] = margin(sys_cl);
    disp("Gm_lqr : ");
    disp(Gm_lqr);
    disp("Pm_lqr : ");
    disp(Pm_lqr);
    disp("Wcg_lqr : ");
    disp(Wcg_lqr);
    disp("Wcp_lqr : ");
    disp(Wcp_lqr);

    % Bode plot
    figure;
    bode(sys_cl);
    title('Bode Plot of LQR-Controlled System');

    % Step response
    figure;
    step(sys_cl);
    title('Step Response of LQR-Controlled System');
    
    
    
end
