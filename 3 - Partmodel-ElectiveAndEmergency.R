## ICU Model Precursor - Part 3
##
## This version uses fixed Poisson-process arrival times
## Random (runif) decision on starting L3 or L2
## Random number of days spent on unit
##
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
  
  
  
  ## Common patient trajectory, where all the timeouts actually occur
  ## runif 0-15 chosen empirically to give a vague balance to the unit
  
  common_patient<-trajectory() %>% 
    timeout(function() {round(runif(1,0,15))}) %>% ## may be L3 or L2 depending on start
    release("halfnurse",function(){ if (get_seized(env,"halfnurse")>1) 1 else 0 }) %>% 
    timeout(function() {round(runif(1,0,15))}) %>% ## definitely L2 or below now
    release_all("halfnurse") %>% 
    release_all("bed")
  
  
  ## Initial trajectory for emergency patients
  ## Appear at timepoint 0.0 so they get priority
  ## If no space, they disappear and count as an emergency transfer out
  
  emergency_transferout<-trajectory() %>% 
    set_global("Emergency Transfer Out",1,mod="+")
  
  emergency_patient<-trajectory() %>% 
    # log_("Emergency Patient") %>% 
    set_attribute("type",1) %>% 
    renege_in(0.7,emergency_transferout) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (runif(1)<0.5) 2 else 1 ## 50% chance of starting at level 3 or 2
    }) %>% 
    renege_abort() %>% 
    join(common_patient)
  
  
  ## Initial trajectory for elective patients
  ## Appear at timepoint 0.1 so they have less priority
  ## If no space, they reappear a week later
  
  elective_delay<-trajectory() %>% 
    set_global("Cancelled Elective Surgery",1,mod="+") %>% 
    timeout(6.5) %>% # 7 days minus the 0.5 before renege
    rollback(3) # to renege_in in elective_patient
  
  elective_patient<-trajectory() %>% 
    #log_("Elective Patient") %>% 
    set_attribute("type",2) %>% 
    timeout(0.1) %>%
    renege_in(0.5,elective_delay) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (runif(1)<0.5) 2 else 1 ## 50% chance of starting at level 3 or 2
    }) %>% 
    renege_abort() %>% 
    join(common_patient)
  
  
  
  ## Display the trajectory for sense-checking (remove before running this multiple times!)  
  #print(plot(elective_patient))
  
  
  ## Set up the simulation. 2 Admission types (for now), 2 resources
  ## bed - represents a physical bed (we have 16)
  ## halfnurse - represents half of a nurse (we have 12 nurses)
  ## halfnurse used because level 2 and below patients can have a 1:2 nursing ratio
  ## nursing coordinators etc should not be included in the numbers
  ## rexp constants empirically derived to be vaguely sensible - 0.33/day for elective, 1/day emergency
  ## rounded so we can use part-days for priority
  
  env %>% 
    add_generator("Emergency Patient",emergency_patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/1)))) %>%  ## Type 1
    add_generator("Elective Patient",elective_patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/3)))) %>%    ## Type 2
    add_resource("bed",16) %>% 
    add_resource("halfnurse",24) %>% 
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
