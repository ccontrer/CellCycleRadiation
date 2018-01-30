using DataStructures

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
  # stab = [sol[:labels][i][curve_point_type]["stab"] for i in 1:len]
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
