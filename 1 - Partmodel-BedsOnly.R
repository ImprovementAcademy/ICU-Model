## ICU Model Precursor - Part 1
##
## This version uses fixed Poisson-process arrival times
## Random (runif) decision on starting L3 or L2
## Random number of days spent on unit
##
## TODO: Add nurses as a resource
## TODO: Split trajectories for elective vs emergency patients
## TODO: Give patients different level-days (use attributes to store requirements). Suggest a common trajectory to generate a patient into the attributes
## TODO: Add generator for patient inter-arrival gaps
## TODO: Change to use real data from CSV
## (c) Tom Lawton, 2019
## Distributed under the GNU General Public License


library("dplyr")
library("parallel")
library("simmer")
library("simmer.plot")
library("ggplot2")
library("lubridate",include.only="round_date") ## Not sure what lubridate is up to but it breaks the sim if the whole thing is included


this.dir <-dirname(parent.frame(2)$ofile)
setwd(this.dir)

## All needs wrapping so that we can run multiple replications (Monte Carlo)

simmer_wrapper <- function(i) {
  
  sim_start_date<-as.Date("2020-01-01")
  sim_end_date<-as.Date("2021-07-01")
  sim_start_num<-as.numeric(sim_start_date)
  sim_end_num<-as.numeric(sim_end_date)
  
  env<-simmer("ICU")
  
  ## Initial trajectory for patients
  ## Appear at timepoint 0.0
  ## If no space, they disappear and count as a transfer out/cancellation
  
  transferout<-trajectory() %>% 
    set_global("Transfer or Cancel",1,mod="+")
  
  patient<-trajectory() %>% 
    renege_in(0.7,transferout) %>% 
    seize("bed") %>% 
    renege_abort() %>% 
    timeout(function() {round(runif(1,0,30))}) %>% ## may be L3 or L2 depending on start
    release_all("bed")
  
  
  ## Display the trajectory for sense-checking (remove before running this multiple times!)  
  #print(plot(patient))
  
  
  ## Set up the simulation. 2 Admission types (for now - using single trajectory), 1 resource
  ## bed - represents a physical bed (we have 16)
  ## rexp constants empirically derived to be vaguely sensible - 0.33/day for elective, 1/day emergency
  ## rounded so we can use part-days for priority
  
  env %>% 
    add_generator("Emergency Patient",patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/1)))) %>%  ## Type 1
    add_generator("Elective Patient",patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/3)))) %>%   ## Type 2
    add_resource("bed",16) %>% 
    run() %>% 
    wrap()
  
}

## Run the simulation once

envs<-simmer_wrapper()

## Some example charts:

## Bed numbers and nursing use for a single run
## nb queueing behaviour is a little odd as a patient can only queue for one resource at a time (until clones are fully implemented)
## therefore whether the queue is for the bed or nurse will depend somewhat on the order they're requested in
print(plot(get_mon_resources(envs),steps=TRUE))
