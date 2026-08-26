using LinearAlgebra
using Statistics
using CSV
using DataFrames
using CIFTI
using GIFTI
using CorticalSurfaces
# using GLMakie
using CairoMakie
using ColorSchemes


networks = Dict(
  0 => nothing, 1 => "DMN", 2 => "VIS", 3 => "FP", 4 => nothing, 5 => "DAN", 
  6=>nothing, 7 => "VAN", 8 => "SAL", 9 => "CO", 10 => "SMD", 11 => "SML", 12 => "AUD", 13 => "Tpole", 14 => "MTL", 15 => "PMN", 16 => "PON")

GordonbrainL = GIFTI.load("brainTemplate/Conte69.L.inflated.32k_fs_LR.surf.gii")
GordonbrainLref = deepcopy(GordonbrainL)
GordonbrainR = GIFTI.load("brainTemplate/Conte69.R.inflated.32k_fs_LR.surf.gii")
Gordonsnp = CIFTI.load("formattedData/roiMedianHeritGordon_GCTAHeatMap.dscalar.nii")[R]
Gordontwin = CIFTI.load("formattedData/roiMedianHeritGordon_twinHeatMap.dscalar.nii")[L]
PFNsnp = CIFTI.load("formattedData/roiMedianHeritPFN_GCTAHeatMap.dscalar.nii")[L]
PFNtwin = CIFTI.load("formattedData/roiMedianHeritPFN_twinHeatMap.dscalar.nii")[L]
θl = -3π / 4  # 45 degrees
rotZr = [cos(θl) -sin(θl) 0;
     sin(θl)  cos(θl) 0;
     0       0      1]
θr = π / 4  # 45 degrees
rotZl = [cos(θr) -sin(θr) 0;
     sin(θr)  cos(θr) 0;
     0       0      1]
GordonbrainL.position.parent.parent .= (rotZl * GordonbrainL.position.parent.parent)
GordonbrainR.position.parent.parent .= (rotZr * GordonbrainR.position.parent.parent)

Label = CIFTI.load("brainTemplate/abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii")
colors = convert.(Int, Label.data)
roiMap = CSV.read("brainTemplate/GRP1_template_parcel.csv", DataFrame)
roilab = map(x -> x==0 ? 0 : roiMap[x, :region], colors)
# convert the colors vector to strings using the networks dictionary 
category_colors = [get(ColorSchemes.rainbow1, i / length(unique(colors))) for i in 1:length(unique(roilab))]  # 4 unique categories
# map from category_colors, but if nothing map to grey rgb value
colorsc = map(x -> x == 0 ? colorant"grey" : category_colors[x], roilab)


function colLimits(vec)
  return(quantile(vec, [0.05, 0.95]))
end



colorMap = :lajolla
def_ax3 = (; :xgridvisible => false,
    :xspinesvisible => false, :yspinesvisible => false, :zspinesvisible => false,
    :xspinesvisible => false, :yspinesvisible => false, :zspinesvisible => false,
    :xgridvisible => false, :ygridvisible => false, :zgridvisible => false,
    :xypanelvisible=> false, :xzpanelvisible => false, :yzpanelvisible=> false,
    :aspect => :data)
def_mesh =(; colormap = colorMap, colorrange = colLimits(Gordonsnp[:,1]), nan_color = :grey, lowclip=:black, highclip = :white)


inch = 96
pt = 4/3
f = Figure(size = (6inch, 6inch), fontsize = 12pt)
Colorbar(f[1,1], limits = colLimits(Gordonsnp[:,1]) * 100, colormap = colorMap, flipaxis = false)
ax = Axis3(f[1, 2], title = "Gordon SNP"; def_ax3...)
mesh!(ax, GordonbrainR, color = Gordonsnp[:,1]; def_mesh...)
hidedecorations!(ax)

ax = Axis3(f[1, 3], title = "Gordon twin"; def_ax3...)
mesh!(ax, GordonbrainL, color = Gordontwin[:,1], colormap = colorMap, colorrange = colLimits(Gordontwin[:,1]))
Colorbar(f[1,4], limits = colLimits(Gordontwin[:,1]) * 100, colormap = colorMap, flipaxis = false)
hidedecorations!(ax)

Colorbar(f[2,1], limits = colLimits(PFNsnp[:,1]) * 100, colormap = colorMap, flipaxis = false)
ax = Axis3(f[2, 2], title = "Personalized SNP"; def_ax3...)
mesh!(ax, GordonbrainR, color = PFNsnp[:,1], colormap = colorMap, colorrange = colLimits(PFNsnp[:,1]))
hidedecorations!(ax)

ax = Axis3(f[2, 3], title = "Personalized twin"; def_ax3...)
hidedecorations!(ax)
mesh!(ax, GordonbrainL, color = PFNtwin[:,1], colormap = colorMap, colorrange = colLimits(PFNtwin[:,1]))
Colorbar(f[2,4], limits = colLimits(PFNtwin[:,1]) * 100, colormap = colorMap, flipaxis = false)

elems = [MarkerElement(color = col, marker = :circle)
  for col in category_colors]
ax = Axis3(f[3, 3], title = "Networks"; def_ax3...)
hidedecorations!(ax)
mesh!(ax, GordonbrainL, color = colorsc, colormap = :rainbow1)
Legend(
  f[3, 4],
  elems,
  [networks[n] for n in unique(colors)];
  markersize = 5, labelsize=8, patchsize = (5,5)
)
# reflect GordonbrainL on x axis
refX = [-1 0 0;
     0  1 0;
     0  0 1]
# GordonbrainLref.position.parent.parent .= (refX * GordonbrainLref.position.parent.parent)
# GordonbrainLref.position.parent.parent .= (rotZr * GordonbrainLref.position.parent.parent)

ax = Axis3(f[3, 2], title = "Networks"; def_ax3...)
hidedecorations!(ax)
mesh!(ax, GordonbrainR, color = colorsc, colormap = :rainbow1)

# Box(f[3, 3])
# Box(f[3, 4])

colgap!(f.layout, 1, Fixed(-150.0))
colgap!(f.layout, 2, Fixed(-250.0))
colgap!(f.layout, 3, Fixed(-150.0))
rowgap!(f.layout, 1, Fixed(-50.0))
rowgap!(f.layout, 2, Fixed(-50.0))
colsize!(f.layout, 1, Fixed(50))
colsize!(f.layout, 2, Fixed(350))
colsize!(f.layout, 3, Fixed(350))
colsize!(f.layout, 4, Fixed(50))
colsize!(f.layout, 4, Fixed(50))
display(f)

save("makieBrain.png", f)

# save colorlab object as a new CIFTI file
colorlab = map(x -> networks[x], roilab)

CIFTI.save("formattedData/networkSurfLabels.dlabel.nii", roilab; template = "brainTemplate/abcd_template_matching_combined_clusters_thresh0.75.dlabel.nii")
