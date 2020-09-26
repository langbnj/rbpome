library(DBI)
library(RMySQL)
superreserveddrv = dbDriver("MySQL")
superreservedcon = dbConnect(superreserveddrv, host="localhost", port=3306, dbname="blang", user="username", pass="password")

Query <- function(query)
{
  a = dbGetQuery(superreservedcon, statement=query)
}
