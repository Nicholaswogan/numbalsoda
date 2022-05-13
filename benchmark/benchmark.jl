
using OrdinaryDiffEq, StaticArrays, BenchmarkTools, Sundials, ModelingToolkit

println("Lorenz")
# Lorenz
function lorenz_static(u,p,t)
    @inbounds begin
        dx = p[1]*(u[2]-u[1])
        dy = u[1]*(p[2]-u[3]) - u[2]
        dz = u[1]*u[2] - p[3]*u[3]
    end
    SA[dx,dy,dz]
end
u0 = SA[1.0,0.0,0.0]
p  = SA[10.0,28.0,8/3]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz_static,u0,tspan,p)
@btime solve(prob,Tsit5(),saveat=0.1,reltol=1.0e-8,abstol=1.0e-8)
@btime solve(prob,Vern8(),saveat=0.1,reltol=1.0e-8,abstol=1.0e-8); 

println("\nRober")
# rober
function rober(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    @inbounds begin
        du[1] = -k₁*y₁+k₃*y₂*y₃
        du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
        du[3] =  k₂*y₂^2
    end
    nothing
end
prob = ODEProblem(rober,[1.0,0.0,0.0],(0.0,1e5),[0.04,3e7,1e4])

@btime solve(prob,Rodas5(),reltol=1.0e-8,abstol=1.0e-8, saveat = 1000)
@btime solve(prob,TRBDF2(),reltol=1.0e-8,abstol=1.0e-8, saveat = 1000) 
@btime solve(prob,CVODE_BDF(),reltol=1.0e-8,abstol=1.0e-8, saveat = 1000) 

# rober
function rober_static(u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du1 = -k₁*y₁+k₃*y₂*y₃
    du2 =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    du3 =  k₂*y₂^2
    SA[du1,du2,du3]
end
prob = ODEProblem{false}(rober_static,SA[1.0,0.0,0.0],(0.0,1e5),SA[0.04,3e7,1e4])
@btime solve(prob,Rodas5(),reltol=1.0e-8,abstol=1.0e-8, saveat = 1000)

function rober_jac(u,p,t)
  y₁,y₂,y₃ = u
  k₁,k₂,k₃ = p
  J11 = k₁ * -1
  J21 = k₁
  J31 = 0
  J12 = y₃ * k₃
  J22 = y₂ * k₂ * -2 + y₃ * k₃ * -1
  J32 = y₂ * 2 * k₂
  J13 = k₃ * y₂
  J23 = k₃ * y₂ * -1
  J33 = 0
  SA[J11 J12 J13
     J21 J22 J23
     J31 J32 J33]
end

ff = ODEFunction(rober_static, jac=rober_jac)
prob2 = ODEProblem{false}(ff,SA[1.0,0.0,0.0],(0.0,1e5),SA[0.04,3e7,1e4])
@btime solve(prob2,Rodas5(), reltol=1.0e-8, abstol=1.0e-8, saveat = 1000); 
