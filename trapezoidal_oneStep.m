function xout = trapezoidal_oneStep(x, dt, A, B)
%  One step of the Trapezoidal method
%
%  Args:
%      x: state vector
%      dt: time step
%      A: state matrix
%      B: independent values vector
%
%  Returns:
%      xout: state vector after one step of the Trapezoidal method


    Aprime = eye(4) - 0.5 * dt * A;
    bprime = (eye(4) + 0.5 * dt * A) * x + dt * B;

    xout  = Aprime\bprime;
    
end