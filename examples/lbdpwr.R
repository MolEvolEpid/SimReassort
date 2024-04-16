simulate("LBDPwr",time=4) |> plot(points=TRUE)

simulate("LBDPwr",lambda=2,mu=1,psi=3,rhoA=.5,rhoB=0,n0=1,time=1) |>
  simulate(time=10,lambda=1) |>
  plot(points=TRUE)

# a pure birth process
simulate("LBDPwr",lambda=1,mu=0,psi=3,rhoA=1,rhoB=0,n0=1,time=4) |>
  plot(points=TRUE)