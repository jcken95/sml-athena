## RUN FIRST ORDER SA FOR EACH MODEL/MODEL COMPONENT

###sml
st.sml <- Sys.time()
saveRDS(st.sml, "start-sml.RDS")
source("sml-mean.R")
source("sml-var.R")
en.sml <- Sys.time()
saveRDS(en.sml, "end-sml.RDS")

###hetgp
st.hgp <- Sys.time()
saveRDS(st.hgp, "start-hgp.RDS")
source("hgp-mean.R")
source("hgp-var.R")
en.hgp <- Sys.time()
print("sml time")
readRDS("end-sml.RDS") - readRDS("start-sml.RDS")
print("hgp time")
en.hgp - readRDS("start-hgp.RDS")