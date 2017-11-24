library(plyr)
library(matrixStats)
library(abind)
options(stringsAsFactors = F)

# Read LiCOR --------------------------------------------------------------


# Function to read Li6400 output files when remarks are mixed with the main data
read_licor_simple = function(path) {
  # Get all the lines of text
  text = readLines(path)
  # Get all indices that start with ", < or $
  comments = grep('^[\\",<,\\$].', text)
  # Lines that do not contain comments
  good_lines = text[-comments]
  # Get the columns with names of variables
  columns = text[grep('^\\$START.',text) + 1]
  # Now we need to separate according to the \t...this is so slow...
  good_lines_split = strsplit(good_lines, "\t")
  # and now we have to coerce it to numeric and put into matrix form...
  output = laply(good_lines_split, function(x) as.numeric(x[-2]))
  # add the names
  names_output = strsplit(columns, '[\\",\t]')[[1]]
  names_output = names_output[-which(names_output == "")]
  names_output = names_output[-2]
  colnames(output) = names_output
  return(output)
}

# Read & process induction curves ----------------------------------------------------------

# Get the names of all the files
raw_data_files = list.files("Input/Experiment/Induction curves/")


get_data <- function(raw_data_files) {
  # genotypes
  genotypes = c("npq1", "col", "npq4","spsa","rca2","rwt43","aba2")
  # biological replicates
  replicates = 1:5
  # the replicates for Col are special because it is the control and was measured at different dates
  replicates_col0 = paste0(rep(1:5, 3),letters[1:3])
  # build the ids for the different measurements
  plants = c(paste(rep(genotypes[-2], each = 5), replicates),
             paste("col",replicates_col0))
  # Create a nested list to store all the data
  # 1st level = Genotype
  # 2nd level = Replicate
  induction_data <- vector("list", length(genotypes))
  names(induction_data) = genotypes
  for(i in names(induction_data)) {
    if(i == "col") {
      use_replicates = replicates_col0
    } else {
      use_replicates = replicates
    }
    induction_data[[i]] = vector("list", length(use_replicates))
    names(induction_data[[i]]) = use_replicates
    for(j in use_replicates) {
      filename = file.path("Input/Experiment/Induction curves/",grep(paste(i,j), raw_data_files, value = TRUE, fixed = TRUE))
      if(length(filename) == 0) {
        induction_data[[i]][[j]] = NA
      } else {
        induction_data[[i]][[j]] = read_licor_simple(filename)
      }
    }
  }
  return(induction_data)
}

induction_data = get_data(raw_data_files)

# Read & process light transients (no darkness) ----------------------------------------------------------

# The measurements corresponding to the light response curves are distributed
# across multiple files

# Each raw data is identified by
# <date> <genotype>-<biological replicate> lrc<part of lrc>
# There may be hyphens or spaces in the file names

# Get the names of all the files
raw_data_files = list.files("Input/Experiment/LRC/")

get_data <- function(raw_data_files) {
  
  # genoypes
  genotypes = c("npq1", "col", "npq4","spsa","rca2","rwt43","aba2")
  replicates = 1:5
  replicates_col0 = paste0(rep(1:5, 3),letters[1:3])
  plants = c(paste(rep(genotypes[-2], each = 5), replicates),
             paste("col",replicates_col0))
  
  # Which are the chunks in which the lrc is divided?
  chunks = c("lrc1","lrc2","lrc3","lrc4")
  
  # Create a nested list to store all the data
  # 1st level = Genotype
  # 2nd level = Replicate
  # 3rd level = Chunk
  
  lrc_data <- vector("list", length(genotypes))
  names(lrc_data) = genotypes
  for(i in names(lrc_data)) {
    if(i == "col") {
      use_replicates = replicates_col0
    } else {
      use_replicates = replicates
    }
    lrc_data[[i]] = vector("list", length(use_replicates))
    names(lrc_data[[i]]) = use_replicates
    for(j in use_replicates) {
      lrc_data[[i]][[j]] = vector("list", 4)
      for(k in 1:4) {
        filename = file.path("Input/Experiment/LRC/",grep(paste(i,j,paste0("lrc",k)), raw_data_files, value = TRUE, fixed = TRUE))
        lrc_data[[i]][[j]][[k]] = read_licor_simple(filename)
      }
    }
  }
  return(lrc_data)
}

lrc_data <- get_data(raw_data_files)

# Fix electronic noise-related errors
lrc_data$col[[9]][[3]][1431,c("Photo","Trmmol","CO2R","H2OR","Cond","CO2S")] = 
  lrc_data$col[[9]][[3]][1430,c("Photo","Trmmol","CO2R","H2OR","Cond","CO2S")]
lrc_data$col[[8]][[4]][1420,c("Photo","Trmmol","CO2R","H2OR","Cond","CO2S")] = 
  lrc_data$col[[8]][[4]][1419,c("Photo","Trmmol","CO2R","H2OR","Cond","CO2S")]

# Name transients
for(i in 1:length(lrc_data)) {
  for(j in 1:length(lrc_data[[i]])) {
    names(lrc_data[[i]][[j]]) = paste0("transient",2:5)
  }
}

# Merge induction and light transients
for(i in 1:length(lrc_data)) {
  for(j in 1:length(lrc_data[[i]])) {
    lrc_data[[i]][[j]][["transient1"]] = induction_data[[i]][[j]]
  }
}
transients = lrc_data
transients = transients
rm(lrc_data, induction_data)


for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      transients[[i]][[j]][[k]] = cbind(transients[[i]][[j]][[k]], Time = transients[[i]][[j]][[k]][,"FTime"])
      transients[[i]][[j]][[k]][,"Time"] = c(1,2,2 + cumsum(diff(transients[[i]][[j]][[k]][-1,"Time"])))
    }
  }
}

# Subset variables required, change names and scales
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      if(k == 5) {
        transients[[i]][[j]][[k]] = transients[[i]][[j]][[k]][,c("Time","Photo","PARi","CO2R","Fm'","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","PhiPS2","VpdL","H2OS","Ci")]
        colnames(transients[[i]][[j]][[k]]) = c("Time","Photo","PAR","CO2R","Fm.","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","PhiPS2","VPD","H2OS","Ci")
      } else {
        transients[[i]][[j]][[k]] = transients[[i]][[j]][[k]][,c("Time","Photo","PARi","CO2R","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","VpdL","H2OS","Ci")]
        colnames(transients[[i]][[j]][[k]]) = c("Time","Photo","PAR","CO2R","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","VPD","H2OS","Ci")
      }
    }
  }
}

# Keep original Cond for later pictures
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      transients[[i]][[j]][[k]] = cbind(transients[[i]][[j]][[k]], oCond = transients[[i]][[j]][[k]][,"Cond"])
    }
  }
}


# Correct negative PAR and very high PAR values
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      transients[[i]][[j]][[k]][,"PAR"] = pmax(transients[[i]][[j]][[k]][,"PAR"], 0)
      weird_PAR = which(transients[[i]][[j]][[k]][,"PAR"] > 1500)
      if(length(weird_PAR) > 0) {
        transients[[i]][[j]][[k]][weird_PAR,"PAR"] = transients[[i]][[j]][[k]][weird_PAR-2,"PAR"]
      }
    }
  }
}

# Correct stomatal conductance (offset after 1 min)
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      if(k == 5) {
        Cond1 = mean(transients[[i]][[j]][[k]][150:170,"Cond"])
        Cond0 = mean(transients[[i]][[j]][[k]][1:20,"Cond"])
        offset = Cond0 - Cond1
        transients[[i]][[j]][[k]][,"Cond"] = c(transients[[i]][[j]][[k]][1:120,"Cond"],
                                               transients[[i]][[j]][[k]][91:120,"Cond"],
                                               transients[[i]][[j]][[k]][-(1:150),"Cond"] + offset)
      } else {
        Cond1 = mean(transients[[i]][[j]][[k]][90:110,"Cond"])
        Cond0 = mean(transients[[i]][[j]][[k]][1:20,"Cond"])  
        offset = Cond0 - Cond1
        transients[[i]][[j]][[k]][,"Cond"] = c(transients[[i]][[j]][[k]][1:60,"Cond"],
                                               transients[[i]][[j]][[k]][31:60,"Cond"],
                                               transients[[i]][[j]][[k]][-(1:90),"Cond"] + offset)     
      }
    }
  }
}

# Smooth measurements (to reduce stiffness of simulations)
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:5) {
      Cond = transients[[i]][[j]][[k]][,c("Time","Cond")]
      Cond_spline = smooth.spline(Cond[,1], Cond[,2],spar = 0.5)
      transients[[i]][[j]][[k]][,"Cond"] = predict(Cond_spline, Cond[,1])$y     
      
      Trmmol = transients[[i]][[j]][[k]][,c("Time","Trmmol")]
      Trmmol_spline = smooth.spline(Trmmol[,1], Trmmol[,2],spar = 0.5)
      transients[[i]][[j]][[k]][,"Trmmol"] = predict(Trmmol_spline, Trmmol[,1])$y
      
      Tleaf = transients[[i]][[j]][[k]][,c("Time","Tleaf")]
      Tleaf_spline = smooth.spline(Tleaf[,1], Tleaf[,2],spar = 0.5)
      transients[[i]][[j]][[k]][,"Tleaf"] = predict(Tleaf_spline, Tleaf[,1])$y   + 273.15      
      
      CO2R = transients[[i]][[j]][[k]][,c("Time","CO2R")]
      CO2R_spline = smooth.spline(CO2R[,1], CO2R[,2],spar = 0.5)
      transients[[i]][[j]][[k]][,"CO2R"] = predict(CO2R_spline, CO2R[,1])$y        
      
      # Detect when PAR changes and after that set PAR constant (i.e. remove small changes)
      par_change = which(abs(diff(transients[[i]][[j]][[k]][,"PAR"])) > 100)
      par_change = par_change[length(par_change)]
      second_light = transients[[i]][[j]][[k]][par_change+2,"PAR"]
      n = nrow(transients[[i]][[j]][[k]])
      transients[[i]][[j]][[k]][(par_change+2):n,"PAR"] = second_light
      
      transients[[i]][[j]][[k]][,"Tair"] = mean(transients[[i]][[j]][[k]][,"Tair"]) + 273.15
    }
  }
}

# Indicate when fluorescence flashes were added
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    k = 5
    flashes = which(abs(diff(transients[[i]][[j]][[k]][,"Fm."])) > 0) + 1
    flashon = rep(0,nrow(transients[[i]][[j]][[k]]))
    flashon[flashes] = 1
    transients[[i]][[j]][[k]] = cbind(transients[[i]][[j]][[k]], flashon = flashon)
    flashes = which(transients[[i]][[j]][[k]][,"flashon"] == 1)
    npq = rep(0, nrow(transients[[i]][[j]][[k]]))
    npq[flashes] = (transients[[i]][[j]][[k]][flashes[1],"Fm."] - transients[[i]][[j]][[k]][flashes,"Fm."])/
      transients[[i]][[j]][[k]][flashes,"Fm."]
    transients[[i]][[j]][[k]] = cbind(transients[[i]][[j]][[k]], NPQ = npq)
  }
}

# Add acclimation
for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in 1:length(transients[[i]][[j]])) {
      new_time = c(0,3600, transients[[i]][[j]][[k]][,"Time"]+3601)
      transients[[i]][[j]][[k]] = rbind(transients[[i]][[j]][[k]][c(1,1),], transients[[i]][[j]][[k]])
      transients[[i]][[j]][[k]][,"Time"] = new_time
    }
  }
}

# Calculate the minimum length for each combination of genotype and type of transient
minmax_len = vector("list",length(transients))
names(minmax_len) = names(transients)
for(i in 1:length(minmax_len)) {
  minmax_len[[i]] = vector("list", 5)
  names(minmax_len[[i]]) = paste0("transient",1:5)
  for(j in 1:length(minmax_len[[i]])) {
    minmax_len[[i]][[j]] = Inf
  }
}

for(i in 1:length(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in names(transients[[i]][[j]])) {
      minmax_len[[i]][[k]] = min(minmax_len[[i]][[k]], nrow(transients[[i]][[j]][[k]]))
    }
  }
}

# Pre-allocate data-frame to store the interpolated transients
transients_all = data.frame(
  Genotype = c(rep("npq1", sum(as.numeric(minmax_len[["npq1"]]))*5),
               rep("col", sum(as.numeric(minmax_len[["col"]]))*15),
               rep("npq4", sum(as.numeric(minmax_len[["npq4"]]))*5),
               rep("spsa", sum(as.numeric(minmax_len[["spsa"]]))*5),
               rep("rca2", sum(as.numeric(minmax_len[["rca2"]]))*5),
               rep("rwt43", sum(as.numeric(minmax_len[["rwt43"]]))*5),
               rep("aba2", sum(as.numeric(minmax_len[["aba2"]]))*5)),
  ID = c(rep(letters[1:5], each = sum(as.numeric(minmax_len[["npq1"]]))),
         rep(letters[1:15], each = sum(as.numeric(minmax_len[["col"]]))),
         rep(letters[1:5], each = sum(as.numeric(minmax_len[["npq4"]]))),
         rep(letters[1:5], each = sum(as.numeric(minmax_len[["spsa"]]))),
         rep(letters[1:5], each = sum(as.numeric(minmax_len[["rca2"]]))),
         rep(letters[1:5], each = sum(as.numeric(minmax_len[["rwt43"]]))),
         rep(letters[1:5], each = sum(as.numeric(minmax_len[["aba2"]])))),
  TransientType = c(rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["npq1"]])),5),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["col"]])),15),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["npq4"]])),5),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["spsa"]])),5),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["rca2"]])),5),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["rwt43"]])),5),
                    rep(rep(paste0("transient",1:5), times = as.numeric(minmax_len[["aba2"]])),5)),
  Obs = NA, Time = NA, Photo = NA, PAR = NA, CO2R = NA, Cond = NA, Trmmol = NA, CO2S = NA, Tleaf = NA,Tair = NA, H2OR = NA,
  Fm. = NA, flashon = NA, NPQ = NA, PhiPS2 = NA, VPD = NA, H2OS = NA, Ci = NA, oCond = NA)

# Store all the data from transients
for(i in names(transients)) {
  for(j in 1:length(transients[[i]])) {
    for(k in names(transients[[i]][[j]])) {
      len = minmax_len[[i]][[k]]
      indices = with(transients_all, which(Genotype == i &
                                             ID == letters[[j]] &
                                             TransientType == k))
      if(k == "transient1")
        transients_all[indices, c("Time","Photo","PAR","CO2R","Fm.","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","flashon","NPQ","PhiPS2","VPD","H2OS","Ci", "oCond")] = 
        transients[[i]][[j]][[k]][1:len,c("Time","Photo","PAR","CO2R","Fm.","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","flashon","NPQ","PhiPS2","VPD","H2OS","Ci", "oCond")]
      else 
        transients_all[indices, c("Time","Photo","PAR","CO2R","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","VPD","H2OS","Ci", "oCond")] = 
        transients[[i]][[j]][[k]][1:len,c("Time","Photo","PAR","CO2R","Cond","Tair","Trmmol","CO2S","Tleaf","H2OR","VPD","H2OS","Ci", "oCond")]
      transients_all[indices, "Obs"] = 1:len
    }
  }
}


# Calculate mean per genotype and transient type
transients_mean = ddply(transients_all, c("Genotype","TransientType","Obs"), function(x) {
  temp = colMeans(x[,-(1:4)])
  if(("flashon" %in% colnames(x)) && !(sum(x[,"flashon"]) %in% c(0,nrow(x)))) temp["NPQ"] = 0
  temp
  })

transients_se = ddply(transients_all, c("Genotype","TransientType","Obs"), 
                      function(x) colSds(as.matrix(x[,-(1:4)]))/sqrt(nrow(x)))

names(transients_se) = names(transients_mean)

# Read & process lightfleck measurements --------------------------------------------

# The measurements corresponding to the light response curves are distributed
# across multiple files

# Each raw data is identified by
# <date> <genotype>-<biological replicate> lf<part of lf>
# There may be hyphens or spaces in the file names

# Get the names of all the files
raw_data_files = list.files("Input/Experiment/Lightflecks/")

get_data <- function(raw_data_files) {
  
  # genotypes
  genotypes = c("npq1", "col", "npq4","spsa","rca2","rwt43","aba2")
  replicates = 1:5
  replicates_col0 = paste0(rep(1:5, 3),letters[1:3])
  plants = c(paste(rep(genotypes[-2], each = 5), replicates),
             paste("col",replicates_col0))
  
  # Which are the chunks in which the lf is divided?
  chunks = c("lf high","lf mid","lf low")
  
  # Create a nested list to store all the data
  # 1st level = Genotype
  # 2nd level = Replicate
  # 3rd level = Chunk
  
  lf_data <- vector("list", length(genotypes))
  names(lf_data) = genotypes
  for(i in names(lf_data)) {
    if(i == "col") {
      use_replicates = replicates_col0
    } else {
      use_replicates = replicates
    }
    lf_data[[i]] = vector("list", length(use_replicates))
    names(lf_data[[i]]) = use_replicates
    for(j in use_replicates) {
      lf_data[[i]][[j]] = vector("list", 3)
      for(k in 1:3) {
        filename = file.path("Input/Experiment/Lightflecks/",grep(paste(i,j,paste("lf",c("low","mid","high")[k])), raw_data_files, value = TRUE, fixed = TRUE))
        lf_data[[i]][[j]][[k]] = read_licor_simple(filename)
      }
    }
  }
  return(lf_data)
}

lf_data <- get_data(raw_data_files)

for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      lf_data[[i]][[j]][[k]] = cbind(lf_data[[i]][[j]][[k]], Time = lf_data[[i]][[j]][[k]][,"FTime"])
      lf_data[[i]][[j]][[k]][,"Time"] = c(1,2,2 + cumsum(diff(lf_data[[i]][[j]][[k]][-1,"Time"])))
    }
  }
}


# Subset variables required, change names and scales

for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      lf_data[[i]][[j]][[k]] = lf_data[[i]][[j]][[k]][,c("Time","Photo","PARi","CO2R","Cond","Tleaf","Trmmol","CO2S","H2OR","Tair")]
      colnames(lf_data[[i]][[j]][[k]]) = c("Time","Photo","PAR","CO2R","Cond","Tleaf","Trmmol","CO2S","H2OR","Tair")
    }
  }
}

# Correct negative PAR and very high PAR values
# Round to nearest 10
for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      lf_data[[i]][[j]][[k]][,"PAR"] = round_any(pmax(lf_data[[i]][[j]][[k]][,"PAR"], 0), 10)
    }
  }
}

# Smooth lightfleck stomatal conductance
for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      Cond = lf_data[[i]][[j]][[k]][,c("Time","Cond")]
      Cond_spline = smooth.spline(Cond[,1], Cond[,2],spar = 0.7)
      lf_data[[i]][[j]][[k]][,"Cond"] = predict(Cond_spline, Cond[,1])$y   
    }
  }
}

# Smooth measurements to remove errors
for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      
      CO2R = lf_data[[i]][[j]][[k]][,c("Time","CO2R")]
      CO2R_spline = smooth.spline(CO2R[,1], CO2R[,2],spar = 0.7)
      lf_data[[i]][[j]][[k]][,"CO2R"] = predict(CO2R_spline, CO2R[,1])$y
      
      lf_data[[i]][[j]][[k]][,"Tair"] = lf_data[[i]][[j]][[k]][,"Tair"] + 273.15
      lf_data[[i]][[j]][[k]][,"Tleaf"] =  lf_data[[i]][[j]][[k]][,"Tleaf"] + 273.15
    }
  }
}

# Add acclimation
for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      
      new_time = c(0,3600, lf_data[[i]][[j]][[k]][,"Time"]+3601)
      lf_data[[i]][[j]][[k]] = rbind(lf_data[[i]][[j]][[k]][c(1,1),], lf_data[[i]][[j]][[k]])
      lf_data[[i]][[j]][[k]][,"Time"] = new_time
      
    }
  }
}

# Create data.frame with average lf_data across genotypes
lf_data_all = data.frame(Genotype = c(rep("npq1", 1922*5*3), 
                                      rep("col", 1922*15*3),
                                      rep(c("npq4","spsa", "rca2", "rwt43", "aba2"), each = 1922*5*3)),
                         ID = c(rep(letters[1:5], each = 1922*3),
                                rep(letters[1:15], each = 1922*3),
                                rep(rep(letters[1:5], each = 1922*3), 5)),
                         Amplitude = c(rep(rep(c(100,200,500),5), each = 1922),
                                       rep(rep(c(100,200,500),15), each = 1922),
                                       rep(rep(rep(c(100,200,500),5), each = 1922), 5)),
                         Obs = rep(1:1922, rep = 135),
                         Time = NA,
                         Photo = NA,
                         PAR = NA,
                         CO2R = NA,
                         Cond = NA,
                         Tleaf = NA,
                         Trmmol = NA,
                         CO2S = NA,
                         H2OR = NA,
                         Tair = NA)

c = 1
for(i in 1:length(lf_data)) {
  for(j in 1:length(lf_data[[i]])) {
    for(k in 1:3) {
      lf_data_all[c:(c + 1921),-(1:4)] = lf_data[[i]][[j]][[k]]
      c = c + 1922
    }
  }
}

lf_data_mean = ddply(lf_data_all, c("Genotype","Amplitude","Obs"), function(x) colMeans(x[,-(1:4)]))

lf_data_se = ddply(lf_data_all, c("Genotype","Amplitude","Obs"), 
                   function(x) colSds(as.matrix(x[,-(1:4)]))/sqrt(nrow(x)))

# Read & process CO2 response curves ------------------------------------------------

# Get the names of all the files
raw_data_files = list.files("Input/Experiment/ACi/")


get_data <- function(raw_data_files) {
  # genotypes
  genotypes = c("npq1", "col", "npq4","spsa","rca2","rwt43","aba2")
  # biological replicates
  replicates = 1:5
  # the replicates for Col are special because it is the control and was measured at different dates
  replicates_col0 = paste0(rep(1:5, 3),letters[1:3])
  # build the ids for the different measurements
  plants = c(paste(rep(genotypes[-2], each = 5), replicates),
             paste("col",replicates_col0))
  # Create a nested list to store all the data
  # 1st level = Genotype
  # 2nd level = Replicate
  aci_data <- vector("list", length(genotypes))
  names(aci_data) = genotypes
  for(i in names(aci_data)) {
    if(i == "col") {
      use_replicates = replicates_col0
    } else {
      use_replicates = replicates
    }
    aci_data[[i]] = vector("list", length(use_replicates))
    names(aci_data[[i]]) = use_replicates
    for(j in use_replicates) {
      filename = file.path("Input/Experiment/ACi/",grep(paste(i,j), raw_data_files, value = TRUE, fixed = TRUE))
      aci_data[[i]][[j]] = read_licor_simple(filename)
    }
  }
  return(aci_data)
}

aci_data <- get_data(raw_data_files)

# Need to append the measurement at high CO2 that repeated separately
aci_aba2_1500 = read_licor_simple("Input/Experiment/ACi/extra_aba2_4_aci")
aci_data[["aba2"]][[4]] = aci_data[["aba2"]][[4]][-c(470:520),]
aci_data[["aba2"]][[4]] = rbind(aci_data[["aba2"]][[4]],
                                aci_aba2_1500[,c(colnames(aci_data[["aba2"]][[4]]))])

## Read the aci leak ##
aci_leak = read.csv("Input/Experiment/aci_leak.csv")

## Round the value of Co2
aci_leak[,1] = c(50, 100, 150, 200, 350, 500, 750, 1000, 1500)

## Sort them in the same order as the aci curve
aci_leak = aci_leak[c(6,5,4,3,2,1,6,7,8,9),]


# Determine the transitions where Fv/Fm changes
# Choose the first point of the transition
get_last_12_points = function(cleaned_data) {
  transition_points = which(abs(diff(cleaned_data)) > 0)
  chunks = vector("list",length(transition_points))
  for(i in 1:length(chunks)) {
    chunks[[i]] = max((transition_points[i]-10),1):min((transition_points[i]+1), length(cleaned_data))
  }
  return(chunks)
}

# Average the last 60 s
vars = c("Genotype","Replicate","CO2", "Photo","Cond","Ci","PhiPS2","PhiCO2","Tleaf","Tair","Trmmol","H2OR")
aci_df <- as.data.frame(matrix(NA, nrow = (length(aci_data) + 2)*5*10, ncol = length(vars)))
names(aci_df) = vars
aci_df[,3] = c(500, 350, 200, 150, 100, 50, 500,750, 1000, 1500)
counter = 0
for(i in 1:length(aci_data)) {
  for(j in 1:length(aci_data[[i]])) {
    # Detect transients based on Fm' changes
    indices_chunks = get_last_12_points(aci_data[[i]][[j]][,"Fm'"])
    # Some manual corrections
    if(i == 2 & j == 15) indices_chunks = indices_chunks[c(1:4,6:11)]
    if(i == 7 & j == 4) indices_chunks = indices_chunks[c(1:9,12)]
    if(i == 7 & j == 3) indices_chunks = indices_chunks[c(6,1:9)]
    # Check that there are no more errors
    if(length(indices_chunks) > 10 | length(indices_chunks) < 10) {
      print(paste0("Error for i = ",i," and j = ",j, collapse = ""))
    }
    # Calculate the means of Photo, Cond, Ci for each CO2 step
    # For LEFII and PhiPS2 simply retrieve the new value
    for(k in 1:10) {
      counter = counter + 1
      aci_df[counter,1:2] = c(names(aci_data)[i], j)
      for(h in c("Photo","Cond","Ci","PhiCO2","Tleaf","Tair","Trmmol","H2OR")) {
        aci_df[counter,h] = mean(aci_data[[i]][[j]][indices_chunks[[k]],h], na.rm = TRUE)
      }
      for(h in c("PhiPS2")) {
        aci_df[counter,h] = aci_data[[i]][[j]][indices_chunks[[k]][length(indices_chunks[[k]])],h]
      }
      # Correct the value of photosynthesis by the leak at this CO2 concentration
      aci_df[counter, "Photo"] = aci_df[counter, "Photo"] - aci_leak[k,2]
    }
  }
}

aci_df[which(aci_df$Genotype == "aba2" & aci_df$Replicate == 3 & aci_df$CO2 == 500)[1],"Photo"] =
  aci_df[which(aci_df$Genotype == "aba2" & aci_df$Replicate == 3 & aci_df$CO2 == 500)[1],"Photo"]*1.0005
# This has created duplications of the values at 500.
# Only the first two are valid, so remove duplications
aci_df = aci_df[!duplicated(aci_df[,1:4]),]


# Average per genotype
aci_df_genotype_avg = ddply(aci_df, c("Genotype", "CO2"),
                            function(x) colMeans(x[,-(1:3)]))

# Standard error per genotype
se_cols = function(data, n) {
  if(n < 3) return(rep(0, length = ncol(data)))
  se = numeric(ncol(data))
  for(i in 1:ncol(data)) {
    se[i] = sd(data[,i])/sqrt(n)
  }
  return(se)
}
aci_df_genotype_se = ddply(aci_df, c("Genotype", "CO2"),
                           function(x) se_cols(x[,-(1:3)],
                                               length(unique(x$Replicate))))

# Put into a single dataframe
colnames(aci_df_genotype_se)[-(1:2)] = paste(colnames(aci_df)[-(1:3)], "se", sep = "_")
aci_df_genotype = cbind(aci_df_genotype_avg, aci_df_genotype_se[,-(1:2)])

# Environmental inputs for A-Ci curve 
aci_inputs = vector("list", length(transients))
names(aci_inputs) = names(transients)
for(i in 1:length(transients)) {
  aci_inputs[[i]] = vector("list", length(transients[[i]]))
  for(j in 1:length(transients[[i]])) {
    subaci = subset(aci_df, Genotype == names(transients)[[i]] & Replicate == j)
    Cond = subaci$Cond
    Ci = subaci$Ci
    CO2 = subaci$CO2
    Tair = subaci$Tair + 273.15 
    Tleaf = subaci$Tleaf + 273.15
    Trmmol = subaci$Trmmol
    H2OR = subaci$H2OR
    aci_inputs[[i]][[j]] = data.frame(Time  = c(0,3600,3601,3600+180,3600+181,3600+180*2,3600+180*2+1,
                                                3600+180*3,3600+180*3+1,3600+180*4,3600+180*4+1,
                                                3600+180*5,3600+180*5+1,3600+180*10,5401,5400+4*60,5400+4*60+1,
                                                5400+8*60,5400+8*60+1,5400+12*60),
                                      Cond  = rep(Cond, each = 2),PAR = rep(1000, 20),Ci = rep(Ci, each = 2),
                                      CO2R = rep(CO2, each = 2), H2OR = rep(H2OR, each = 2), Tair = rep(Tair, each = 2),
                                      Tleaf = rep(Tleaf, each = 2), Trmmol = rep(Trmmol, each = 2),
                                      Genotype = rep(names(transients)[[i]], 20),
                                      Replicate = rep(j,20), Obs = 1:20)
  }
}

# Aggregate by genotype
for(i in 1:length(aci_inputs)) {
  for(j in 1:length(aci_inputs[[i]])) {
    if(i == 1 & j == 1)
      aci_inputs_all = aci_inputs[[i]][[j]]
    else
      aci_inputs_all = rbind(aci_inputs_all,aci_inputs[[i]][[j]])
  }
}

aci_inputs_means = ddply(aci_inputs_all, c("Genotype","Obs"), 
                         function(x) colMeans(x[,c("Time","Cond","PAR","Ci","CO2R","H2OR","Tleaf","Tair","Trmmol")]))

aci_inputs_se = ddply(aci_inputs_all, c("Genotype","Obs"), 
                      function(x) colSds(as.matrix(x[,c("Time","Cond","PAR","Ci","CO2R", "H2OR", "Tleaf","Tair","Trmmol")]))/sqrt(nrow(x)))
names(aci_inputs_se) = names(aci_inputs_means)

# Save all processed files and empty memory -------------------------------

save(aci_df_genotype, aci_inputs_means, aci_inputs_se , aci_df, file = "Intermediate/aci.RData")
save(lf_data_mean, lf_data_se, file = "Intermediate/lightflecks.RData")
save(transients_mean, transients_se, transients_all, file = "Intermediate/transients.RData")


rm(list = ls())
gc()
