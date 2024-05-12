
function error = emgfittingObjective(params, x, Y)
    al = params(1);
    mu = params(2);
    x0 = params(3);
    s = params(4);
    b = params(5);

    model = (al/2) * exp((s^2 ./ (2 * x0^2)) - ((x - mu)/x0)) .* erfc (((s^2) - (x0*(x-mu))) / (2^(1/2) * s *x0)) + b;

    error = sum((model - Y).^2);
end
