module Rocket3D
# using PyPlot
using Plots
# Plots.plotly()
# Plots.plotlyjs()
Plots.pyplot()


const G = 6.67408e-11
const rEarth = 12742000/2 #m
const mEarth = 5.972e24 #kg
const vEarthRotation = 460 #m/s

type FlightPlan
    fireTimes::Array{Float64, 2}
    thrustDirection::Array{Float64, 2}
end
type Body
    emptyMass::Float64
    fuelMass::Float64
    thrust::Float64
    fuelPerSecond::Float64
    flightPlan::FlightPlan

    r::Vector{Float64}
end

function getMass(rocket::Body)
    return rocket.emptyMass + rocket.fuelMass
end

function main()

    Plots.close("all")

    #I will define the flight plan in terms of spherical coordinates
    flightPlan = FlightPlan(zeros(3,2), zeros(3,3))
    part=1 #go straight up
    flightPlan.fireTimes[part,1] = 0.0
    flightPlan.fireTimes[part,2] = 300.0
    flightPlan.thrustDirection[part,1] = 1.0 #r
    flightPlan.thrustDirection[part,2] = 0.0 #theta
    flightPlan.thrustDirection[part,3] = 0.0 #phi

    part=2 #45 degree thrust
    flightPlan.fireTimes[part,1] = 300.0
    flightPlan.fireTimes[part,2] = 700.0
    flightPlan.thrustDirection[part,1] = 1.0/sqrt(2.0) #r
    flightPlan.thrustDirection[part,2] = 1.0/sqrt(2.0) #theta
    flightPlan.thrustDirection[part,3] = 0.0 #phi

    part=3 #only in theta direction to gain orbital speed
    flightPlan.fireTimes[part,1] = 1000.0
    flightPlan.fireTimes[part,2] = 1040.0
    flightPlan.thrustDirection[part,1] = 0.0 #r
    flightPlan.thrustDirection[part,2] = 1.0 #theta
    flightPlan.thrustDirection[part,3] = 0.0 #phi

    startPosition = [rEarth,0.0, 0.0]
    rocket = Body(24.0, 120.0, 1600.0, 0.1, flightPlan, startPosition)

    timestep= 1.0
    time, x, y, z = runCalc(timestep, rocket)

    plt1 = Plots.plot(time, getR(x, y, z) - rEarth,
                    markershape = :hexagon,
                    markersize = 2,
                    markeralpha = 0.6,
                    markercolor = :blue,
                    markerstrokecolor = :blue,
                    label="Timestep = $(timestep)",
                    xlabel = "Time [s]",
                    ylabel = "Height [m]",
                    reuse=false)
    Plots.display(plt1)
    Plots.savefig("TimeVsHeight_dt=$(timestep).jpg")

    plt2 = Plots.plot(y, x,
                    markershape = :hexagon,
                    markersize = 2,
                    markeralpha = 0.6,
                    markercolor = :blue,
                    markerstrokecolor = :blue,
                    label="Timestep = $(timestep)",
                    xlabel = "y [m]",
                    ylabel = "x [m]",
                    reuse = false)

    # timestep= 0.01
    # rocket.fuelMass = 120.0
    # rocket.r = startPosition
    # time, x, y, z = runCalc(timestep, rocket)

    # Plots.plot!(plt2, y, x,
    #             markershape = :circle,
    #             markersize = 2,
    #             markercolor = :yellow,
    #             markerstrokecolor = :yellow,
    #             label="Earth surface",
    #             label="Timestep = $(timestep)")

    # plt = path3d(1, #xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
    #             xlab = "x", ylab = "y", zlab = "z",
    #             #title = "Lorenz Attractor",
    #             marker = 1,
    #             reuse = false)
    # for i=1:size(x)[1]
    #   push!(plt, x[i], y[i], z[i])
    # end
    # Plots.display(plt)


    # plt3 = Plots.plot(earthY, earthX)
    # Plots.display(plt3)



    earthX, earthY = getEarth()
    Plots.plot!(plt2, earthY, earthX,
                markershape = :hexagon,
                markersize = 2,
                markercolor = :red,
                markerstrokecolor = :red,
                label="Earth surface",
                xlim = (minimum(y)-10^3, maximum(y)+10^3),
                ylim = (minimum(x)-10^3, maximum(x)+10^3))

    Plots.display(plt2)
    Plots.savefig("X-Y_timestep=$(timestep).jpg")
end

function getEarth()
    earthX=[]
    earthY=[]
    for θ1 = 0:0.001:2.0*π
        append!(earthX, rEarth*sin(θ1))
        append!(earthY, rEarth*cos(θ1))
    end
    return earthX, earthY
end

function getXYZfromR(r, θ, φ)
    x = r.*cos.(θ).*sin.(φ)
    y = r.*sin.(θ).*sin.(φ)
    z = r.*cos.(φ)
    return x, y, z
end

function getR(r)
  return sqrt.(r[1].^2 + r[2].^2 + r[3].^2)
end

function getR(x, y, z)
  return sqrt.(x.^2 + y.^2 + z.^2)
end
function getRθφFromXYZ(x, y, z)
  r = sqrt.(x.^2 + y.^2 + z.^2)
  θ = atan.(y./x)
  φ = acos.(z./r)
  return r, θ, φ
end

function getRDirection(r)
    r, θ, φ = getRθφFromXYZ(r[1], r[2], r[3])
    return getRDirection(r, θ, φ)
end
function getRDirection(r, θ, φ)
    r, θ, φ = [cos.(θ)*sin.(φ), sin.(θ)*sin.(φ), cos.(φ)]
end

#r is x, y, z, returns theta direction
function getθDirection(r)
    r, θ, φ = getRθφFromXYZ(r[1], r[2], r[3])
    return getθDirection(r, θ, φ)
end
function getθDirection(r, θ, φ)
    return [-sin(θ), cos(θ), 0]
end

#r is x, y, z, returns phi direction
function getφDirection(r)
    r, θ, φ = getRθφFromXYZ(r[1], r[2], r[3])
    return getφDirection(r, θ, φ)
end
function getφDirection(r, θ, φ)
    return [cos(θ)*cos(φ), sin(θ)*cos(φ), -sin(φ)]
end

function calculateThurstDirection(time::Float64, r::Array{Float64, 1}, fp::FlightPlan)
    #Check ascent profile
    dir = [0, 0, 0]
    ft= fp.fireTimes
    for i = 1:size(ft)[1]
      if ft[i, 1] <= time && ft[i, 2] > time

        r, θ, φ = getRθφFromXYZ(r[1], r[2], r[3])
        dirR = getRDirection(r, θ, φ)
        dirθ = getθDirection(r, θ, φ)
        dirφ = getφDirection(r, θ, φ)
        # println("r, θ, φ = $(dirR), $(dirθ), $(dirφ) ")
        # println("fp.thrustDirection[i, 1] = $(fp.thrustDirection[i, :]) ")
        dir =  dirR*fp.thrustDirection[i, 1] +
                dirθ*fp.thrustDirection[i, 2] +
                dirφ*fp.thrustDirection[i, 3]
        # println("dir = $(dir) ")
        break
      end
    end
    return dir
end

function getGravitationalForce(r, m1, m2)
    magr = getR(r)
    f = -G*m1*m2/magr^2
    dirR = getRDirection(r)
    return f*dirR
end

function runCalc(δt::Float64, rocket::Body)
    t=0.0
    v=0.0

    xHist = Float64[]
    yHist = Float64[]
    zHist = Float64[]
    time = Float64[]

    c = 0
    println("Starting loop ... ")

    append!(time, t)
    append!(xHist, rocket.r[1] )
    append!(yHist, rocket.r[2])
    append!(zHist, rocket.r[3])

    tic()
    while( getR(rocket.r) >= rEarth-1.0 && t<5000)

        c += 1
        force = vcat(0.0, 0.0, 0.0) #Working in r, θ, φ basis now
        if rocket.fuelMass > 0.0
            dir = calculateThurstDirection(t, rocket.r, rocket.flightPlan)

            force = rocket.thrust * dir
            rocket.fuelMass = rocket.fuelMass - rocket.fuelPerSecond*δt
        end
        # print("force = $(force), ")
        forceG = getGravitationalForce(rocket.r, mEarth, getMass(rocket))
        force =  force + forceG

        if abs(getR(rocket.r) - rEarth < 1.0 && getR(force) < 0.0) #If we do not have enough thrust to take off, stay grounded
            force = vcat(0.0, 0.0, 0.0)
        end

        if c%10 == 0
            print("c = $(c) (t=$(t) s): h=$(getR(rocket.r) - rEarth), r=$(rocket.r)")
            print("force = $(force), fuel = $(rocket.fuelMass)")
            print("\n")
        end
        a = force/getMass(rocket)
        v = v + a*δt

        rocket.r = rocket.r + v*δt
        t = t + δt

        append!(time, t)
        append!(xHist, rocket.r[1] )
        append!(yHist, rocket.r[2])
        append!(zHist, rocket.r[3])

        # println(" force = $(force) mass = $(mass) massFule = $(fuelMass) ")
    end
    println("\nFinished loop after $(c) steps, flight time $(c*δt) s in $(toq()) s of CPU time.")
    return time, xHist, yHist, zHist
end

main()
end
