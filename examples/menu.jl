using REPL.TerminalMenus

examples_dir = joinpath(@__DIR__)
example_files = filter(f -> splitext(f)[1] != "menu", sort(readdir(examples_dir)))
example_names = [splitext(f)[1] for f in example_files]

menu = RadioMenu(example_names, pagesize = 10)

println("Select an example to run:")
choice = request(menu)

if choice == -1
    println("No selection made. Exiting.")
else
    file = joinpath(examples_dir, example_files[choice])
    rel_path=joinpath("examples", example_files[choice])
    println("\nRunning: $rel_path")
    println(repeat("─", length(rel_path)+9))
    include(file)
end