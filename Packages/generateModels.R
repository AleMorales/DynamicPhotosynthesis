library(SimulationModels)

# Original
compile_and_install(ode_file = "MiniModel.ode",name_model = "MiniModel", 
                    output = "MiniModel",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# Original - Ci as input
compile_and_install(ode_file = "MiniModelCi.ode",name_model = "MiniModelCi", 
                    output = "MiniModelCi",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# Emulate rwt43
compile_and_install(ode_file = "MiniModelrwt43.ode",name_model = "MiniModelrwt43", 
                    output = "MiniModelrwt43",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")


# Emulate rwt43 - Ci as input
compile_and_install(ode_file = "MiniModelrwt43.ode",name_model = "MiniModelrwt43", 
                    output = "MiniModelrwt43",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# Emulate rwt43 - Ci as input
compile_and_install(ode_file = "MiniModelrwt43Ci.ode",name_model = "MiniModelrwt43Ci", 
                    output = "MiniModelrwt43Ci",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS qE
compile_and_install(ode_file = "MiniModelQssqE.ode",name_model = "MiniModelQssqE", 
                    output = "MiniModelQssqE",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS qM
compile_and_install(ode_file = "MiniModelQssqM.ode",name_model = "MiniModelQssqM", 
                    output = "MiniModelQssqM",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

Sys.setenv(HOME="C:/Users/Alejandro")
# QSS qI
compile_and_install(ode_file = "MiniModelQssqI.ode",name_model = "MiniModelQssqI", 
                    output = "MiniModelQssqI",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS fR
compile_and_install(ode_file = "MiniModelQssfR.ode",name_model = "MiniModelQssfR", 
                    output = "MiniModelQssfR",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS fRB
compile_and_install(ode_file = "MiniModelQssfRB.ode",name_model = "MiniModelQssfRB", 
                    output = "MiniModelQssfRB",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS PR
compile_and_install(ode_file = "MiniModelQssPR.ode",name_model = "MiniModelQssPR", 
                    output = "MiniModelQssPR",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")

# QSS qs
compile_and_install(ode_file = "MiniModelQssgs.ode",name_model = "MiniModelQssgs", 
                    output = "MiniModelQssgs",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")


# QSS - All QSS in the above combined
compile_and_install(ode_file = "MiniModelQss.ode",name_model = "MiniModelQss", 
                    output = "MiniModelQss",unit_analysis = "true", install = FALSE, 
                    directory = "Code/Models")


# Test the model ----------------------------------------------------------

# Induction curve
library(MiniModelQss)
model = generate_MiniModelQss_model()
model$set_settings("atol", 1e-10)
model$set_time(1:(8e3*3))
sim = cvode(model)


model$set_forcings("Ir", cbind(c(1,1000,1001,2000,2001,3000,3001,4000,4001,5000,5001,6000,6001,7000,7001,8000)*3, 
                               c(0,0, 50, 50, 100, 100, 300, 300, 500, 500, 1000, 1000, 1500, 1500, 2000, 2000)))
model$set_forcings("Ib", cbind(c(1,2), c(0,0)))
model$set_forcings("Ci", cbind(c(1,1000), c(400,400)))
load("Intermediate/ParametersLiterature.RData")
load("Intermediate/ParametersExperiment.RData")
parameters = params_literature
for(i in row.names(params_experiment)) {
  parameters[i,] = params_experiment[i,]
}
filter = which(row.names(parameters) %in% names(model$Parameters$Values))
model$set_parameters(row.names(parameters[filter,]),unname(parameters[filter,"value"]))

test = cvode(model)

plot(test[(1:8)*1e3*3,c("PhiII","NPQ")], t = "o")
plot(test[(1:8)*1e3*3,c("qP","PhiIIo")], t = "o")
plot(test[(1:8)*1e3*3,c("PhiIIo","NPQ")], t = "o")
plot(test[(1:8)*1e3*3,c("qP","PhiqE")], t = "o")
