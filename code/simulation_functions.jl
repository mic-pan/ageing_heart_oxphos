"""
Run multiple simulations of model with different parameters. The parameters should be given as an array of pairs (e.g. [[a => 1, b => 1],[a => 1, b => 2]])
"""
function multi_sim(sys,tspan,ps;u0=[],reltol=1e-7,abstol=1e-10,alg=Rodas5())
    sols = []
    for p in ps
        prob = ODEProblem(sys,u0,tspan,p)
        push!(sols,solve(prob,alg; reltol=reltol, abstol=abstol))
    end
    return sols
end

"""
Extract variables from a solution. The 
"""
extract_var(sols,p) = [sol.retcode==:Success ? sol[p][end] : missing for sol in sols]

