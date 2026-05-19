# SPDX-FileCopyrightText: 2022 Uwe Fechner
# SPDX-License-Identifier: MIT

# Create a Julia system image including all project packages.
# Run via: bin/create_sys_image

using Pkg

# Ensure PackageCompiler is available in the base environment
if !haskey(Pkg.project().dependencies, "PackageCompiler")
    Pkg.add("PackageCompiler")
end

using PackageCompiler

project_dir = dirname(@__DIR__)
sysimage_path = joinpath(project_dir, "sys_MTKCircuits.so")
precompile_script = joinpath(project_dir, "examples", "coupled_system.jl")

packages = [
    :ControlPlots,
    :LaTeXStrings,
    :ModelingToolkit,
    :ModelingToolkitStandardLibrary,
    :OrdinaryDiffEq,
    :OrdinaryDiffEqFIRK,
]

println("Creating system image at: $sysimage_path")
println("Packages: $(join(string.(packages), ", "))")

create_sysimage(
    packages;
    sysimage_path = sysimage_path,
    precompile_execution_file = precompile_script,
    project = project_dir,
)

println("System image created: $sysimage_path")
println("Use it with: julia --sysimage $sysimage_path --project")
