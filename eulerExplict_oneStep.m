function xout = eulerExplict_oneStep(func,t, x, dt)
%   One step of the Euler Explicit method
%
%   Args:
%       func: function to be integrated
%       t: time
%       x: state vector
%       dt: time step
%
%   Returns:
%       xout: state vector after one step of the Euler Explicit method

    xout = x + dt * func(t,x);

end