function [x] = conjgrad(A,b)

TOL = 0.000001;
ITERATION = 1000;

iniX = zeros(size(b));

residual = b - A * iniX;

jacobi = eye(size(A))*4;

h = jacobi\residual;


% direction = -h;

% delta = residual' * h; - 
delta = norm(residual) / norm(b);

it = 0;
while(delta > TOL && it < ITERATION)
%       disp(h);
    a = A*h;
    temp = a' * h;
    
    lambda = (residual' * h)/(h'*A*h);
    
    
    
    nextX = iniX + lambda * h;
    
    % update x and residual;
    iniX = nextX;
    residual = residual - lambda * a;
    
    %update P
    temp2 = jacobi \ residual;
    h = temp2 - ((temp2'*a) / temp) * h;

    
    delta = norm(residual) / norm(b);
    
    
    it = it + 1;
end

disp(it);
x = iniX;

end