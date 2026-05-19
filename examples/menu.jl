using REPL.TerminalMenus

examples_dir = joinpath(@__DIR__)
example_files = filter(f -> splitext(f)[1] != "menu" && f != "LocalPreferences.toml", sort(readdir(examples_dir)))
example_names = [[splitext(f)[1] for f in example_files]; "quit"]

while true
    local menu = RadioMenu(example_names, pagesize = 10)

    println("Select an example to run (or 'quit' to exit):")
    local choice = request(menu)

    if choice == -1 || choice == length(example_names)
        println("Exiting.")
        break
    else
        local file = joinpath(examples_dir, example_files[choice])
        local rel_path = joinpath("examples", example_files[choice])
        println("\nRunning: $rel_path")
        println(repeat("─", length(rel_path) + 9))
        include(file) |> display
        println()
    end
end