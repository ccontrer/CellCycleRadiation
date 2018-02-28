using DataStructures, NamedTuples

type BifurcationCurve{DType,SPType}
  d::DType
  stab::Vector{String}
  special_points::SPType
  changes::Vector{Int}
end

function bifurcation_curve_ccc(PC,name)  # Get the curve
  sol = PC[:curves][name][:sol]
  pts = PyDict(sol)
  d = OrderedDict{Symbol,Vector{Float64}}()
  for k in keys(pts)
    d[Symbol(k)] = pts[k]
  end
  len = length(d[first(keys(d))])

  # Get the stability
  # S => Stable
  # U => Unstable
  # N => Neutral

  # Get this from the information at:
  # https://github.com/robclewley/pydstool/blob/master/PyDSTool/PyCont/ContClass.py#L218
  #TODO: unify stability for all kind of curves
  # curve_point_type = "EP"
  #if calc_stab
  stab = String[]
  for i in 1:len
      p_found = false
      for (key, val) in sol[:labels][i]
          p_stab = get(val, "stab", nothing)
          if p_stab != nothing
              if !p_found
                  push!(stab, p_stab)
                  p_found = true
              else
                  warn("Multiple stability labels found")
              end
          end
      end
  end
  #  stab = [sol[:labels][i][curve_point_type]["stab"] for i in 1:len]
  #else
  #  stab = []
  #end

  curve = PC[:curves][name]
  # Get information for special points, ex limit points
  special_points = Dict{String,Any}()
  for k in keys(curve[:BifPoints])
    special_points[k] = []
    for i in 1:length(curve[:BifPoints][k][:found])
      tmp_dict = PyDict(curve[:BifPoints][k][:found][i]["X"])
      dd = Dict{Symbol,Float64}()
      for k2 in keys(tmp_dict)
        dd[Symbol(k2)]=tmp_dict[k2]
      end
      push!(special_points[k], dd)
    end
  end

  changes = find_changes(stab)

  BifurcationCurve(d,stab,special_points,changes)
end

function find_changes(stab)
  changes = Int[]
  for i in 2:length(stab)
    if stab[i]!= stab[i-1]
      push!(changes,i)
    end
  end
  changes
end

function ntupdate(x::NamedTuple, y::Dict)
    z = @NT()
    for key in keys(x)
        if key in keys(y)
            z = NamedTuples.setindex(z, key, y[key])
        else
            z = NamedTuples.setindex(z, key, x[key])
        end
    end
    z
end
               
function find_limit_cycle(f, u0, tspan, params; alg=Rosenbrock23(), kwargs...)
    function division(u,t,integrator)
        ## Parameters
        ## Set of condition: MPF crosses threshold from positive to negative
        u[1] - params[31]#params.theta_M
    end # ModelCond
    function cycle_stop!(integrator)
        terminate!(integrator)
    end
    cb_start = ContinuousCallback(division, nothing, cycle_stop!)
    # Integrate until first division
    # TODO: assert tspan is larger than needed
    prob = ODEProblem(f, u0, tspan, params)
    sol = solve(prob, alg, callback=cb_start; kwargs...)
    # One full cycle
    u0_start = sol.u[end]
    prob = ODEProblem(f, u0_start, tspan, params)
    sol = solve(prob, alg, callback=cb_start; kwargs...)
    # Done
    sol
end #function find_cell_cycle

function follow_limit_cycle(initcycle, pdomain, stepsize; kwargs...)
    #p[32] should be p.Mass
    cycle = deepcopy(initcycle)
    f = cycle.prob.f
    p = [cycle.prob.p[key] for key in keys(cycle.prob.p)]
    # Start
    minmax = extrema(cycle[1, :])
    T = cycle.t[end] - cycle.t[1]
    out = [[p[32], minmax[1], minmax[2], T]]
    # Forward
    while (p[32] < (pdomain[2] - stepsize))
        #TODO: Fix this, use copy or something
        p[32] += stepsize
        u0 = cycle.u[end]
        tspan = (0.0, 1.3*T)
        cycle = find_limit_cycle(f, u0, tspan, p; kwargs...)
        minmax = extrema(cycle[1, :])
        T = cycle.t[end] - cycle.t[1]
        push!(out, [p[32], minmax[1], minmax[2], T])
    end
    # Backward
    cycle = deepcopy(initcycle)
    f = cycle.prob.f
    p = [cycle.prob.p[key] for key in keys(cycle.prob.p)]
    while ((pdomain[1] + stepsize) < p[32])
        p[32] -= stepsize
        u0 = cycle.u[end]
        tspan = (0.0, 1.8*T)
        cycle = find_limit_cycle(f, u0, tspan, p; kwargs...)
        minmax = extrema(cycle[1, :])
        T = cycle.t[end] - cycle.t[1]
        push!(out, [p[32], minmax[1], minmax[2], T])
    end
    temp = sort(out, by = x -> x[1])
    temp = hcat(temp...)
    dict1 = Dict([:min=>temp[2,:], :max=>temp[3,:], :length=>temp[4]])
    dict2 = OrderedDict(Dict([:Mass=>temp[1,:], :MPF=>dict1]))
    BifurcationCurve(dict2, repeat(["LC"], inner=[length(out)]), Dict{String,Any}(), Int[])
end
