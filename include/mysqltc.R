library(RMariaDB)

if (exists("superreservedcon"))
{
  try(dbDisconnect(superreservedcon), silent=T)
  rm(superreservedcon)
  rm(superreserveddrv)
}

superreserveddrv = dbDriver("MariaDB")

superreservedcon = dbConnect(superreserveddrv, host="localhost", port=3306, dbname="blang", user="username", pass="password", bigint="numeric")

Query <- function(query)
{
  a = dbGetQuery(superreservedcon, statement=query)
}
