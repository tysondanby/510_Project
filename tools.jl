
function newtonitter(x,f)
    dx = 1e-7
    dfdx = (f(x+dx)-f(x-dx))/(2*dx)
    return x - f(x)/dfdx
end

function zeronewton(f,xi)
    xnew = copy(xi)
    xold = 2*xnew + 1
    while abs((xnew-xold)/xold) > 1e-3
        #println(xold)
        #println(xnew)
        xold = copy(xnew)
        xnew = newtonitter(xnew,f)
    end
    return xnew
end