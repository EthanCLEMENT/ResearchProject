function TF_with_feedback
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
   
    % ntf = tf([2000],[260 200  2943])


    % Transfer function computation
    sys = ss(A, B, C, D);
    transferFunction = tf(sys);
    
    %bode(transferFunction);
    [Gm, Pm, Wcg, Wcp] = margin(sys);
    disp("Gm : ");
    disp(Gm);
    disp("Pm : ");
    disp(Pm);
    disp("Wcg : ");
    disp(Wcg);
    disp("Wcp : ");
    disp(Wcp);

    % PID Controller Parameters
    Kp = 200;  % Proportional Gain
    Ki = 50;  % Integral Gain
    Kd = 20;  % Derivative Gain

    pidController = pid(Kp,Ki,Kd);

    % Closed-loop transfer function with PID controller
    sys_cl_pid = feedback(sys*pidController,1);

    [Gm_pid, Pm_pid, Wcg_pid, Wcp_pid] = margin(sys_cl_pid);
    disp("Gm_pid : ");
    disp(Gm_pid);
    disp("Pm_pid : ");
    disp(Pm_pid);
    disp("Wcg_pid : ");
    disp(Wcg_pid);
    disp("Wcp_pid : ");
    disp(Wcp_pid);
    
    figure
    bode(sys);
    figure
    bode(sys_cl_pid);
    figure
    subplot(2,1,1);
    impulse(transferFunction, '-');
    title('Impulse response transfer function');
    
    subplot(2,1,2);
    step(transferFunction, '-');
    title('Step response transfer function');

    figure
    subplot(2,1,1);
    step(sys_cl_pid, '-');
    title('Closed loop PID step ');

    subplot(2,1,2);
    impulse(sys_cl_pid, '-');
    title('Closed loop PID impulse');
   
end
