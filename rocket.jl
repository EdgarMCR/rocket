# using Plots
# pyplot() # Switch to using the PyPlot.jl backend
using PyPlot


G = 6.67408e-11
rEarth = 12742000/2
mEarth = 5.972e24

TIMESTEP = 1

function getGravitationalForce(r, m1, m2)
    -G*m1*m2/r^2
end

function runCalc(δt)
    emptyMass = 30
    fuelMass = 140
    mass = emptyMass + fuelMass
    thrust = 1000
    fuelPerSecond = 5
    r = 12742000/2
    v = 0
    t=0

    pos = zeros(0)
    time = zeros(0)
    c =0
    println("Starting loop ... ")
    append!(pos, r-rEarth)
    append!(time, t)
    while(r >= rEarth-1)
        c += 1
        # if c%100 ==0
        #     print(" step $(c), r= $(r), v = $(v)")
        # end
        force = 0
        if fuelMass > 0
            force = thrust
            fuelMass = fuelMass - fuelPerSecond*δt
            mass = mass - fuelPerSecond*δt
        end

        force += getGravitationalForce(r, mEarth, mass)
        if abs(r - rEarth) < 1 && force < 0 #If we do not have enough thrust to take off, stay grounded
            force = 0
        end

        a = force/mass
        v = v + a*δt
        r = r + v*δt
        t = t + δt
        append!(pos, r-rEarth)
        append!(time, t)
        # println(" force = $(force) mass = $(mass) massFule = $(fuelMass) ")
    end
    println("\nFinished loop after $(c) steps.")
    return time, pos
end

PyPlot.close("all")
fig, ax = PyPlot.subplots(figsize=(10,8))

TIMESTEP = 1.0
time, pos = runCalc(TIMESTEP)
ax[:plot](time, pos, marker = true, markeredgecolor="k", markersize=3, lw=1,label = "$(TIMESTEP)")
# ax[:plot](time, pos , linewidth=1, marker = 'o', markersize=2, label = "dt = $(TIMESTEP)")

TIMESTEP = 0.01
time, pos = runCalc(TIMESTEP)
ax[:plot](time, pos, marker = true, markeredgecolor="k", markersize=3, lw=1, label = "$(TIMESTEP)")
# ax[:plot](time, pos , linewidth=1, marker = 'o', markersize=2,label = TIMESTEP)
# ax[:ylabel]="Height [m]"
# ax[:xlabel]="Time [s]"
# ax[:legend](loc="best")
# display(h)
PyPlot.xlabel("Time [s]")
PyPlot.ylabel("Height [m]")

# PyPlot.title("Your Title Goes Here")
PyPlot.legend(loc="best")
PyPlot.grid("on")
savefig( "figure1")
