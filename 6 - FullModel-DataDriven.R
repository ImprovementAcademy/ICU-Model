## ICU Model v1.0 - Part 6 of teaching
##
## Now uses real data to drive model
## Totals calculated to drive poisson statistic for generator
## Patient pulled from data to drive days at different levels for trajectories
## Uses mclapply so we can run Monte-Carlo
##
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


## Read in patient data

patients<-read.csv("ICNARC Test Data.csv",header=TRUE,sep=",",colClasses=c("Admission.Date"="character"))
patients$Admission.Date<-as.Date(patients$Admission.Date,format="%d/%m/%Y")


## Calculate totals for later use in Poisson distribution

patients$weekday<-factor(weekdays(patients$Admission.Date))
totals<-patients %>% group_by(Admission.Type,weekday) %>% summarise(count=n())
daterange<-difftime(max(patients$Admission.Date),min(patients$Admission.Date),units="days")


## All needs wrapping so that we can run multiple replications (Monte Carlo)

simmer_wrapper <- function(i) {
  
  sim_start_date<-as.Date("2020-01-01")
  sim_end_date<-as.Date("2021-07-01")
  
  env<-simmer("ICU")
  
  
  ## Trajectory to pull a real patient from the list
  ## Presently just chooses one of the same type of patient
  
  choose_real_patient<-trajectory() %>% 
    set_attribute(c("lvl_start","l3days","l2days","l1days","l0days"),function() {
      type<-switch(get_attribute(env,"type"),"Unplanned local admission","Planned local surgical admission","Repatriation","Unplanned transfer in")
      patient<-patients[sample(which(patients$Admission.Type==type),1),]
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
    join(choose_real_patient) %>% 
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
    join(choose_real_patient) %>% 
    timeout(0.1) %>% 
    renege_in(0.5,elective_delay) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (get_attribute(env,"lvl_start")==3) 2 else 1
    }) %>% 
    renege_abort() %>% 
    join(common_patient)
  
  
  ## Initial trajectory for repatriation patients
  ## Appear at timepoint 0.2 so they have less priority
  ## If no space, they reappear a day later
  
  repatriation_rejected<-trajectory() %>% 
    set_global("Repatriation Rejected",1,mod="+") %>% 
    timeout(0.7) %>% # try again tomorrow = 1 day minus 0.3 before renege
    rollback(3) # to renege_in in repatriation_patient
  
  repatriation_patient<-trajectory() %>% 
    #log_("Repatriation Patient") %>% 
    set_attribute("type",3) %>% 
    join(choose_real_patient) %>% 
    timeout(0.2) %>% 
    renege_in(0.3,repatriation_rejected) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (get_attribute(env,"lvl_start")==3) 2 else 1
    }) %>% 
    renege_abort() %>% 
    join(common_patient)
  
  
  ## Initial trajectory for transferred-in patients
  ## Appear at timepoint 0.3 so they have lowest priority
  ## If no space, they are rejected and disappear
  
  transferin_rejected<-trajectory() %>% 
    set_global("Transfer In Rejected",1,mod="+")
  
  transferin_patient<-trajectory() %>% 
    #log_("Transfer In Patient") %>% 
    set_attribute("type",4) %>% 
    join(choose_real_patient) %>% 
    timeout(0.3) %>% 
    renege_in(0.2,transferin_rejected) %>% 
    seize("bed") %>% 
    seize("halfnurse", function() {
      if (get_attribute(env,"lvl_start")==3) 2 else 1
    }) %>% 
    renege_abort() %>% 
    join(common_patient)
  
  
  ## Function to generate time-gaps between patients
  ## Uses closures so it can be replicated across different patient types, and to persist variables across calls
  ## Produces patients from sim_start_date to sim_end_date before returning -1 to signal the end
  ## Patient numbers (per day) are simply a Poisson distribution indexed by day of the week and admission type 
  
  patient_gen<-function(type) {
    date<-sim_start_date
    patients_today<-0
    last_patient<-0
    function() {
      if (patients_today>0){
        patients_today<<-patients_today-1
        return(0)
      } else {
        repeat {
          patients_today<<-rpois(1,max(totals$count[(totals$weekday==weekdays(date))&(totals$Admission.Type==type)]*7/as.numeric(daterange),0))
          if (patients_today>0){
            patients_today<<-patients_today-1
            gap<-date-last_patient
            last_patient<<-date
            date<<-date+1
            return(gap)
          } else {
            date<<-date+1
            if(date>=sim_end_date){return(-1)}
          }
        }
      }
    }
  }
  
  
  ## Set up the simulation. 4 Admission types, 2 resources
  ## bed - represents a physical bed (we have 16)
  ## halfnurse - represents half of a nurse (we have 12 nurses)
  ## halfnurse used because level 2 and below patients can have a 1:2 nursing ratio
  ## nursing coordinators etc should not be included in the numbers
  
  env %>% 
    add_generator("Emergency Patient",emergency_patient,patient_gen("Unplanned local admission")) %>%         ## Type 1
    add_generator("Elective Patient",elective_patient,patient_gen("Planned local surgical admission")) %>%    ## Type 2
    add_generator("Repatriation Patient",repatriation_patient,patient_gen("Repatriation")) %>%                ## Type 3
    add_generator("Transfer In Patient",transferin_patient,patient_gen("Unplanned transfer in")) %>%          ## Type 4
    add_resource("bed",16) %>% 
    add_resource("halfnurse",24) %>% 
    run() %>% 
    wrap()
  
}


runs<-1
#runs<-24

## Run the simulation repeatedly (Monte-Carlo method)
envs<-mclapply(1:runs,simmer_wrapper)

## Some example charts:

## Bed numbers and nursing use for a single run
## nb queueing behaviour is a little odd as a patient can only queue for one resource at a time (until clones are fully implemented)
## therefore whether the queue is for the bed or nurse will depend somewhat on the order they're requested in
print(plot(get_mon_resources(envs[1]),steps=TRUE))

## What is our distribution of markers like?
## Use max values as simmer records each time the value changes, we're only interested in the final value
attribs<-get_mon_attributes(envs)
max_attribs<-attribs %>% group_by(key,replication) %>% summarise(value=max(value))
print(ggplot(max_attribs,aes(x=key,y=value,fill=key)) +
        geom_boxplot() +
        stat_summary(fun.y=mean, geom="point", shape=23,size=6) +
        theme_bw(base_size=16))


## What is our distribution of bed use like?
## nb - system includes queue, server doesn't
resources<-get_mon_resources(envs)
resources2<-dplyr::filter(resources,resource=="bed")
resources2$date<-as.Date(resources2$time,origin="1970-01-01")
resources2$rwdate<-round_date(resources2$date,unit="week")
print(ggplot(resources2,aes(x=date,y=system,color=replication)) +
        geom_point(alpha=0.1,shape=16) +
        scale_color_gradient(low="blue", high="red") +
        stat_summary(aes(x=rwdate,y=system),fun.data="mean_sdl",geom="smooth",se=TRUE) +
        labs(x="Date",y="Beds") +
        theme_bw(base_size=16))


