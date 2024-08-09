##' Congo Ebola Incidence data from the 2018 Ebola epidemic in the Democratic
##' Republic of Congo.
##'
##' @format ## `Congo.EbolaIncidence`
##' \describe{
##'   \item{SitReptNo.}{An integer labelling the report from the reporting
##'                     agency.}
##'   \item{Date}{The date of the report in YYYY-MM-DD format.}
##' }
##'
##' Further columns are named for the location in which the integer number of
##' infections ocurred.
##'
##' \preformatted{
##'   # A tibble: 75 × 30
##'      SitReptNo. Date                Alimbongo  Beni Biena Butembo  Goma Kalunguta Katwa Kayna Kyondo Lubero Mabalako
##'           <dbl> <dttm>                  <dbl> <dbl> <dbl>   <dbl> <dbl>     <dbl> <dbl> <dbl>  <dbl>  <dbl>    <dbl>
##'    1          1 2018-08-05 00:00:00         0     3     0       2     0         0     0     0      0      0       34
##'    2          2 2018-08-12 00:00:00         0     2     0       0     0         0     0     0      0      0       11
##'    3          3 2018-08-20 00:00:00         0     1     0       0     0         0     0     0      0      0       37
##'    4          4 2018-08-26 00:00:00         0     5     0       0     0         0     0     0      0      0        3
##'    5          5 2018-09-02 00:00:00         0     8     0       0     0         1     0     0      0      0        1
##'    6          6 2018-09-09 00:00:00         0     5     0       2     0         0     0     0      0      0        1
##'    7          7 2018-09-16 00:00:00         0     5     0       3     0         0     0     0      0      0        2
##'    8          8 2018-09-23 00:00:00         0     4     0       1     0         0     0     0      0      0        1
##'    9          9 2018-10-02 00:00:00         0    10     0       1     0         0     0     0      0      0        0
##'   10         10 2018-10-07 00:00:00         0    14     0       4     0         0     0     0      0      0        1
##'   # ℹ 65 more rows
##'   # ℹ 17 more variables: Manguredjipa <dbl>, Masereka <dbl>, Musienene <dbl>, Mutwanga <dbl>, Nyiragongo <dbl>,
##'   #   Oicha <dbl>, Pinga <dbl>, Vuhovi <dbl>, Ariwara <dbl>, Bunia <dbl>, Komanda <dbl>, Lolwa <dbl>, Mambasa <dbl>,
##'   #   Mandima <dbl>, Nyakunde <dbl>, Rwampara <dbl>, Tchomia <dbl>
##'   # ℹ Use `print(n = ...)` to see more rows
##' }
"Congo.EbolaIncidence"

##' Reporting health zones of the Democratic Republic of Congo
##'
##' @format ## `healthZonesCongo`
##' \describe{
##'   \item{HealthZone}{The name of the Congolese health zone}
##'   \item{Latitude}{The geographical latitude of the health zone}
##'   \item{Longitude}{The geographical longitude of the health zone}
##' }
##'
##' \preformatted{
##'        HealthZone   Latitude Longitude
##'   1     Alimbongo -0.3655150  29.19118
##'   2          Beni  0.4911300  29.47306
##'   3         Biena  0.5792300  29.11563
##'   4       Butembo  0.1406920  29.33501
##'   5          Goma -1.6582710  29.22013
##' }
"healthZonesCongo"

##' 2018 Congo Ebola epidmeic data
##'
##' "Seed" data to use for a simulation of the 2018 Ebola epidemic in the Congo
##'
##' Here is an example the data used.
##'
##' \preformatted{
##'   Location Latitude Longitude Vaccinated Exposed Infected Recovered Dead
##'   Beni     0.49113  29.47306  0          24      12       0         4
##'   Butembo  0.140692 29.335014 0          0       0        0         0
##'   Mabalako 0.461257 29.210687 0          0       0        0         0
##'   Mandima  1.35551  29.08173  0          0       0        0         0
##' }
##'
##' @format ## `initialInfections.fourCities`
##' \describe{
##'   \item{Location}{The name of the location (reporting health zone) which the other variables describe}
##'   \item{Latitude}{The geographical latitude of the health zone}
##'   \item{Longitude}{The geographical longitude of the health zone}
##'   \item{InitialVaccinated}{The number of vaccinated people}
##'   \item{InitialExposed}{The number of exposed people}
##'   \item{InitialInfections}{The number of infected people}
##'   \item{InitialRecovered}{The number of recovered people}
##'   \item{InitialDead}{The number of dead people}
##' }
"initialInfections.fourCities"
