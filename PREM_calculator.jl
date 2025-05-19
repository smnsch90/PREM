# calculating density as function of Radius using PREM (Dziewonski & Anderson, 1981)
# smn | 07.01.2025
# assuming spherical symmetry and neglecting centrifugal effects

using Plots

# define density function from DA81
function density(r)
       x = r / 6371.0  # normalized grid
       ρ = zeros(length(r))  # initialize density array
       for i = 1:length(r)
              if r[i] < 1221.50  # inner core
                     ρ[i] = 13.0885 - 8.831 * x[i]^2
              elseif r[i] < 3480.00  # outer core
                     ρ[i] = 12.5815 - 1.2638 * x[i] - 3.6426 * x[i]^2 - 5.5281 * x[i]^3
              elseif r[i] < 5701.00  # lower mantle
                     ρ[i] = 7.95650 - 6.4761 * x[i] + 5.5283 * x[i]^2 - 3.0807 * x[i]^3
              elseif r[i] < 5771.00  # transition zone 1
                     ρ[i] = 5.3197 - 1.4836 * x[i]
              elseif r[i] < 5971.00  # transition zone 2
                     ρ[i] = 11.2494 - 8.0298 * x[i]
              elseif r[i] < 6151.00  # transition zone 3
                     ρ[i] = 7.1089 - 3.8045 * x[i]
              elseif r[i] < 6346.60  # LVZ + LID
                     ρ[i] = 2.6910 + 0.6924 * x[i]
              elseif r[i] < 6356.00  # lower crust
                     ρ[i] = 2.900
              elseif r[i] < 6368.00  # upper crust
                     ρ[i] = 2.600
              else  # ocean
                     ρ[i] = 1.020
              end
       end

       return ρ * 1e3  # convert density to kg/m^3
end

# define mass function
function mass(r, ρ)
       avgρ = zeros(length(r))     # initialize average density array
       M = zeros(length(r))        # initialize mass array
       r = r * 1e3                 # convert radius to meters
       for i = 2:length(r)
              avgρ[i] = (ρ[i] + ρ[i-1]) / 2                                  # average density - trapezoidal rule
              M[i]    = M[i-1] + 4/3 * π * (r[i]^3 - r[i-1]^3) * avgρ[i]     # mass
       end

       return M
end

# define gravity function
function gravity(r, M)
       G = 6.67430e-11      # Gravitational constant in m^3 kg^-1 s^-2
       g = zeros(length(r)) # initialize gravity array
       r = r * 1e3          # convert radius to meters
       for i = 2:length(r)
              g[i] = G * M[i] / r[i]^2  # gravity
       end

       return g
end

# define pressure function
function pressure(r, ρ, g)
       revρ    = reverse(ρ)        # reverse density array
       revg    = reverse(g)        # reverse gravity array
       avgrevρ = zeros(length(r))  # initialize average density array
       avgrevg = zeros(length(r))  # initialize average gravity array
       P = zeros(length(r))        # initialize pressure array
       r = r * 1e3                 # convert radius to meters
       for i = 2:length(r)
              avgrevρ[i] = (revρ[i] + revρ[i-1]) / 2           # average density - trapezoidal rule
              avgrevg[i] = (revg[i] + revg[i-1]) / 2           # average gravity - trapezoidal rule
              P[i]       = P[i-1]   + avgrevρ[i] * avgrevg[i] * (r[i] - r[i-1])   # pressure
       end

       return P * 1e-9  # convert pressure to GPa
end

# constants
R = 6371.0  # Earth's radius in km

# numerics
n  = 1e4        # number of grid points
dr = R / (n-1)  # grid spacing in km
r  = 0:dr:R     # grid in km

ρ = density(r)       # call density function
M = mass(r, ρ)       # call mass function
g = gravity(r, M)    # call gravity function
P = pressure(r, ρ, g)# call pressure function

# plot results as subplots
p1 = plot(reverse(r), ρ/1e3, xlabel="Depth [km]", ylabel="Density [kg/m³ ⋅ 10³]", title="PREM density profile", legend=false, linewidth=2, linecolor=:black)
p2 = plot(reverse(r), M / 1e24, xlabel="Depth [km]", ylabel="Mass [kg ⋅ 10^{24}]", title="PREM mass profile", legend=false, linewidth=2, linecolor=:black)
p3 = plot(reverse(r), g, xlabel="Depth [km]", ylabel="Gravity [m/s²]", title="PREM gravity profile", legend=false, linewidth=2, linecolor=:black)
p4 = plot(r, P, xlabel="Depth [km]", ylabel="Pressure [GPa]", title="PREM pressure profile", legend=false, linewidth=2, linecolor=:black)

plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 800), 
       titlefont=font(14, "Arial"), 
       guidefont=font(12, "Arial"), 
       tickfont=font(10, "Arial"), 
       framestyle=:box,
       margins = 25Plots.px,
       )
       
# export figure as png
#savefig("PREM_densities.png")

# export as pdf
#savefig("PREM_densities.pdf")