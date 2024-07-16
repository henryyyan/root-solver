function [root, info] = modifiedzeroin3038415876(func, Int, params)
    % Check if the initial interval has a sign change
    info = struct;
    info.flag = 0;
    root = NaN;
    a = Int.a;
    b = Int.b;
    f = func;
    array = []; % keep previous values of f(x3)
    root_tol = params.root_tol;
    func_tol = params.func_tol;
    maxit = params.maxit;

    if abs(func(a)) <= func_tol
        root = a;
        return;
    end

    if abs(func(b)) <= func_tol
        root = b;
        return;
    end

    if func(a) * func(b) >= 0
        info.flag = 1;
        return;
    end

    x0 = a;
    x1 = b;
    x2 = (a + b) / 2;

    consecutive = 0;

    for iter = 1:maxit
        if b - a <= root_tol % if root_tol met
            root = (b + a) / 2;
            return;
        end
        % Inverse Quadratic Interpolation
        x3 = IQI(f, x0, x1, x2);
       
        % if x3 not in [a, b]
        if x3 < a || x3 > b
            % perform bisection
            [anew, bnew] = bisection(f, a, b);
           
            if abs(f(anew)) <= func_tol
                root = anew;
                return
            end
            if abs(f(bnew)) <= func_tol
                root = bnew;
                return
            end
            
            % check if new interval has sign change
            if f(anew) * f(bnew) < 0
                a = anew;
                b = bnew;
                x0 = a;
                x1 = b;
                x2 = (a + b) / 2;
               
                % Reset consecutive counter
                consecutive = 0;
            else
                info.flag = 1
                return;
            end
          
        else % Using IQI
            if abs(f(x3)) <= func_tol
                root = x3;
                return;
            end
            consecutive = consecutive + 1;
            array(end + 1) = f(x3);

            if consecutive >= 4
                previous = array(end - 3); % consider consecutive = 4
                if abs(previous / f(x3)) < 2
                    % Perform Bisection steps
                    [anew, bnew] = bisection(f, a, b);
                    a = anew;
                    b = bnew;
                    x0 = a;
                    x1 = b;
                    x2 = (a + b) / 2;
                    
                    consecutive = 0;
                else
                    x0 = x1;
                    x1 = x2;
                    x2 = x3;
                end
        
            else
                % Update values for the next iteration
                x0 = x1;
                x1 = x2;
                x2 = x3;
            end
        end
    end

    % Maximum number of iterations reached
    root.flag = 1;
    return;
end

function x3 = IQI(f, x0, x1, x2)
    f0 = f(x0);
    f1 = f(x1);
    f2 = f(x2);

    x3 = (x0 * f1 * f2) / ((f0 - f1) * (f0 - f2)) + (x1 * f0 * f2) / ((f1 - f0) * (f1 - f2)) + (x2 * f0 * f1) / ((f2 - f0) * (f2 - f1));
end

function [anew, bnew] = bisection(f, a, b)
    c = (a + b) / 2;
    f_new = f(c);
   
    if f(a) * f_new < 0
        anew = a;
        bnew = c;
    else
        anew = c;
        bnew = b;
    end
end
