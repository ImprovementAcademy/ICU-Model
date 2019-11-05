## ICU Model Precursor - Part 4
##
## This version uses fixed Poisson-process arrival times
## Now generates proper L3/2/1/0 days - stored using attributes
##
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




## All needs wrapping so that we can run multiple replications (Monte Carlo)

simmer_wrapper <- function(i) {
  
  sim_start_date<-as.Date("2020-01-01")
  sim_end_date<-as.Date("2021-07-01")
  sim_start_num<-as.numeric(sim_start_date)
  sim_end_num<-as.numeric(sim_end_date)
  
  env<-simmer("ICU")
  
  
  ## Trajectory to create a fake patient
  ## Random number of days of each type
  ## Has to store the starting level so we know how many halfnurses to seize
  
  create_fake_patient<-trajectory() %>% 
    set_attribute(c("lvl_start","l3days","l2days","l1days","l0days"),function() {
      patient<-data.frame(
        Level.3.days=rpois(1,5),
        Level.2.days=rpois(1,5),
        Level.1.days=rpois(1,3),
        Level.0.days=rpois(1,1)
      )
      lvl_start<-0
      if (patient$Level.1.days>0) lvl_start<-1
      if (patient$Level.2.days>0) lvl_start<-2
      if (patient$Level.3.days>0) lvl_start<-3
      return(c(lvl_start,patient$Level.3.days,patient$Level.2.days,patient$Level.1.days,patient$Level.0.days))
    })
  
  
  ## Common patient trajectory, where all the timeouts actually occur
  ## Initially there's an offset to account for the way ICNARC dates are calculated (eg stay 8pm-8am, counts as 2)
  ## Offset subtraction means that patients will "leave" on the morning of their last day to let the bed be re-used
  
  common_patient<-trajectory() %>% 
    set_attribute(c("l0days","l1days","l2days","l3days"),function(){
      offset<-now(env)-floor(now(env))+0.9
      ldays<-c(get_attribute(env,"l0days"),get_attribute(env,"l1days"),get_attribute(env,"l2days"),get_attribute(env,"l3days"))
      for (i in 1:length(ldays)) {
        red<-min(ldays[i],offset)
        ldays[i]<-ldays[i]-red
        offset<-offset-red
      }
      return(ldays)
    }) %>%
    timeout(function() {get_attribute(env,"l3days")}) %>% 
    release("halfnurse",function(){ if (get_seized(env,"halfnurse")>1) 1 else 0 }) %>% 
    timeout(function() {get_attribute(env,"l2days")}) %>%
    ## Could put delayed discharge logic in here
    timeout(function() {get_attribute(env,"l1days")}) %>%
    timeout(function() {get_attribute(env,"l0days")}) %>% 
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
    join(create_fake_patient) %>% 
    renege_in(0.7,emergency_transferout) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (get_attribute(env,"lvl_start")==3) 2 else 1
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
    join(create_fake_patient) %>% 
    timeout(0.1) %>%
    renege_in(0.5,elective_delay) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (get_attribute(env,"lvl_start")==3) 2 else 1
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
    add_generator("Emergency Patient",emergency_patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/1)))) %>%
    add_generator("Elective Patient",elective_patient,from_to(sim_start_num,sim_end_num,function() round(rexp(1,1/3)))) %>%
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
