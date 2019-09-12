Valid =[
    5.5384    0.7026  119.8596;
    21.6027   12.8981  108.9038;
    36.8732   23.1884   95.1366;
    39.8395   16.3060   92.0798;
    25.4536    6.7102  111.5937;
];
% fun = @SquareErr;
fun = @Err;

options = optimoptions('fsolve','Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1e10,'MaxIterations',1e10,'StepTolerance',1e-10,'FunctionTolerance',1e-10);

[x,fval,exitflag,output] = fsolve(fun,Valid,options);