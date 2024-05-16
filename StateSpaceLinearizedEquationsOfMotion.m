function computeAndPlotOutput
    % Parameters
    c = 0.5;   % Damping coefficient
    m = 0.8;   % Mass
    g = 9.81;  % Gravitational acceleration
    dc = 1;    % Distance from pivot to center of mass
    d1 = 2;    % Distance from pivot to thrust point
    m_motor = 0.2; % Mass of motor
    J = 1/3*(m*d1^2) + m_motor*d1^2; % Moment of inertia

    % Define the linearized system matrices A and B
    A = [0, 1; (-m*g*dc)/J, -c/J];
    B = [0; d1/J];

    % Output matrix C
    C = [1, 0];

    % D matrix
    D = 0;

    % Set initial conditions in degrees
    initial_theta = 0;         % Initial angle (in degrees)
    initial_theta_dot = 0;     % Initial angular velocity

    initial_conditions = [deg2rad(initial_theta); deg2rad(initial_theta_dot)];

    % Time span for simulation
    tspan = [0 45];

    % Input function (thrust torque)
    u = @(t) 1 * d1;

    % Simulate the linearized system
    [t, x] = ode45(@(t, x) linearPendulumODE(t, x, A, B, u), tspan, initial_conditions);
    x(:,1) = rad2deg(x(:,1));

    % Compute the output y = Cx + Du
    y = (C * x' + D * u(t)').';

    % Display eigenvalues of A
    disp("Eigenvalues of A:");
    disp(eig(A));

    % Check if the controllability matrix has full rank
    disp("Check if rank Ctrb(A,B) has full rank: ");
    CtrbMatrice = ctrb(A,B);
    disp(rank(CtrbMatrice));

    % Plot the output of the linearized system
    figure;
    subplot(2,1,1);
    plot(t, y, 'LineWidth', 1);
    title('Output of Linearized System (Output: Angle)');
    xlabel('Time (s)');
    ylabel('Output (Angle)');

    % LQR Control Design
    Q = [10, 0; 0, 1];  % State cost weights
    R = 1;              % Control cost weight
    K = lqr(A, B, Q, R);% Compute LQR gain matrix

    % Display LQR gain matrix
    disp("LQR Gain K: ");
    disp(K);

    % Define reference state (30 degrees, 0 angular velocity)
    ref_theta = deg2rad(30);    % Reference angle (in radians)
    ref_theta_dot = 0;           % Reference angular velocity

    % State feedback control law with reference state
    Ac = A - B * K;
    Bc = B;
    Cc = C;
    Dc = D;

    % Simulate the LQR-controlled system with reference state
    [t_lqr, x_lqr] = ode45(@(t, x) LQRlinearPendulumODE(t, x, Ac, Bc, u, ref_theta), tspan, initial_conditions);
    x_lqr(:,1) = rad2deg(x_lqr(:,1));

    % Compute the output for the LQR-controlled system
    y_lqr = (Cc * x_lqr' + Dc * u(t_lqr)').';

    % Plot the LQR-controlled output
    subplot(2,1,2);
    plot(t_lqr, y_lqr, 'LineWidth', 1);
    title('Output of LQR-controlled System (Output: Angle)');
    xlabel('Time (s)');
    ylabel('Output (Angle)');
end

function dxdt = linearPendulumODE(t, x, A, B, u)
    % Linearized ODE function for the pendulum with thrust
    dxdt = A * x + B * u(t);
end

function dxdt = LQRlinearPendulumODE(t, x, Ac, Bc, u, ref_theta)
    % Linearized ODE function for the pendulum with LQR control
    % Compute control input based on the error between current state and reference state
    error = ref_theta - x(1);  % Error in angle
    u_control = -Ac(2, :) * x + Bc * u(t) - Ac(2, 1) * ref_theta; % Control law
    dxdt = Ac * x + Bc * u_control;
end
