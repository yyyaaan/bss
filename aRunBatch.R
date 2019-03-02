library(bigrquery) 
library(parallel)
library(magrittr)
reset_access_cred()
set_service_token("/home/yanpan/.gcp.json")
bqcon <- dbConnect(bigrquery::bigquery(), project = "yyyaaannn", dataset = "BSS", billing = "yyyaaannn")

paste0(getwd(), "/sim/")



df_test <- data.frame(x = 1, y = 99)


tryCatch(dbWriteTable(bqcon, "test", df_test), 
         error = function(e) print(paste(e, "it goes!")))

parallel_save <- function(){
  mclapply(10:50,
           function(vec) sim_bootstrap(100*2^{0:10}, vec, savefile = "res_PAR_"),
           mc.cores = detectCores() - 1)
}

