
# Function to create a package that contains the C++ implementation of the model and the R code
# that generates an ODEmodel or DAEmodel
create_package = function(R_files, Cpp_files, name, path = ".", overwrite = TRUE, license = "MIT") {

  project_path = file.path(path, name)

  if(overwrite) unlink(project_path, recursive = TRUE)

  Rcpp.package.skeleton(name = name,path = path,force = overwrite, cpp_files = file.path(path, Cpp_files),
                        code_files = R_files, attributes = TRUE, environment = environment(),
                        example_code = FALSE)

  description = read.dcf(file.path(project_path, "DESCRIPTION"))
  description[,"LinkingTo"] = "RcppSundials, RcppArmadillo, Rcpp"
  description <- cbind(description, Imports = "RcppSundials, RcppArmadillo, Rcpp")
  description <- cbind(description, Depends = "SimulationModels")
  write.dcf(x = description, file = file.path(project_path, "DESCRIPTION"))

  file.create(file.path(project_path, "src/Makevars"))
  write("CXX_STD = CXX11", file = file.path(project_path, "src/Makevars"))

}

# Wrapper on the command to install the package in the system
install_package = function(name, path = ".") {
  project_path = file.path(path, name)
  print(paste0("R CMD INSTALL ", project_path, collapse = " "))
  system(paste0("R CMD INSTALL ", project_path, " --no-multiarch", collapse = " "))
}

# Create and install the package (combines the two functions in the above)
create_and_install_package = function(model_inputs, Cpp_files, name, solver = "RcppSundials", path = ".",
                                      overwrite = TRUE, author = "Your Name", email = "your@email.com", license = "MIT") {
  R_files = create_R_files(model_inputs, name, path, solver)
  create_package(R_files, Cpp_files, name, path, overwrite, author, email, license)
  install_package(name, path)
}


# Compile the ode files in Julia. Create and install the package (combines the two functions in the above)
compile_and_install = function(ode_file = "MiniModelChamber.ode", name_model = "MiniModelChamber",
                              output = "MiniModelChamber",directory = getwd(),
                              unit_analysis = "true",solver = "RcppSundials", overwrite = TRUE,
                              remove = TRUE, install = TRUE) {
  #curdir = getwd()
  #setwd(directory)
  compile_ode(ode_file = ode_file, name_model = name_model,
                       output = output, directory = directory,
                       unit_analysis = unit_analysis)
  model_inputs = paste0(output,"_inputs.json",collapse = "")
  Cpp_files = paste0(output,".cpp",collapse = "")
  R_files = create_R_files(model_inputs, name_model, directory, solver)
  create_package(R_files, Cpp_files, name_model, directory, overwrite)
  if(remove) {
     file.remove(paste0(file.path(directory, output),".cpp",collapse = ""))
     file.remove(paste0(file.path(directory, output),".R",collapse = ""))
     file.remove(paste0(file.path(directory, output),"_inputs.json",collapse = ""))
  }
  if(install) install_package(name_model, directory)
  #setwd(curdir)
}

# Create the R files for the model
# It assumes that the derivatives are available and that the names
# were generated automatically by ODEDSL.
# Also, it needs to specify the library to perform integration (RcppSundials or deSolve)
create_R_files = function(model_inputs, name, directory, solver = "RcppSundials") {

  # Load the inputs from the JSON file
  inputs = fromJSON(file = file.path(directory, model_inputs))

  # With this information, the code generation can be constructed.
  # I use string interpolation with @{<code>} and use the create input function below that relies on paste
code = qq('generate_@{name}_model = function() {
  ODEmodel$new(
    @{create_input("States",inputs$state, inputs$coef_states)},
    @{if("forcings" %in% names(inputs))
          paste0(create_input("Forcings",inputs$forcings, inputs$coef_forcings),",",collapse = "")}
    @{if("parameters" %in% names(inputs))
          paste0(create_input("Parameters",inputs$parameters, inputs$coef_parameters),",",collapse = "")}
    @{if("observed" %in% names(inputs))
          paste0(create_input("Observed",inputs$observed, inputs$coef_observed),",",collapse = "")}
    model_function = \'getNativeSymbolInfo(name = "@{name}",PACKAGE = "@{name}")$address\',
    Settings = list(rtol = 1e-6, atol = 1e-6, method = "bdf",positive = 1)
  )
}')

  # Create the R file and write the code into it
  f = file(file.path(directory, qq('@{name}.R')),"w")
  cat(code, file = f)
  close(f)
  return(file.path(directory, qq('@{name}.R')))
}

# Construct the code necessary to generate an input as Input = list(Value = c(name = value), Coefs = c(name = value))
create_input = function(name, values, coefs) {
  if(name == "Observed") {
  qq('@{name} = list(Names = c(@{paste0(\'"\',values,\'"\' ,collapse = ", ")}),
   Coefs = c(@{paste(names(coefs), " = ", coefs, collapse = ", ")}))')
  } else if (name == "Forcings") {
    forctimes = vector("list",length(values))
    forcvalues = vector("list",length(values))
    if(length(forcvalues) > 0) {
      for(i in 1:length(forcvalues)) {
        forctimes[[i]] = values[[i]][[1]]
        forcvalues[[i]] = values[[i]][[2]]
      }
    }
  qq('@{name} = list(Values = list(@{paste0(names(values), " = cbind(", forctimes, ", " , forcvalues, collapse = "), ")})),
   Coefs = c(@{paste(names(coefs), " = ", coefs, collapse = ", ")}))')
  } else {
  qq('@{name} = list(Values = c(@{paste0(names(values), " = ", values, collapse = ", ")}),
   Coefs = c(@{paste(names(coefs), " = ", coefs, collapse = ", ")}))')
  }
}


# Compile the code by opening Julia and running basic script
compile_ode = function(ode_file = "MiniModelChamber.ode", name_model = "MiniModelChamber",
                       output = "MiniModelChamber", directory = getwd(),
                       unit_analysis = "true") {
  #curdir = getwd()
  #setwd(directory)
  language = "cpp"
  julia_code = qq('julia -e \'cd(\"@{directory}\"); import legacyODEDSL; ode = legacyODEDSL; source = \"@{ode_file}\"; ode_model =  ode.OdeModel(ode.OdeSource(ode.process_file(source))); ode.generate_code!(ode_model, language = \"@{language}\", unit_analysis = @{unit_analysis}, name = \"@{name_model}\", file = \"@{output}\",jacobian = false);\'')
  julia_code =  gsub('"',"\\\"",julia_code, fixed = TRUE)
  julia_code =  gsub("'",'"', julia_code, fixed = TRUE)
  cat(julia_code,'\n')
  suppressWarnings({
  message = system(julia_code)
  })
  if(message == 1) stop("The ODE file could not be processed correctly. Please see output above for details.")
  #setwd(curdir)
}







