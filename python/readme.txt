These are notes from 12/31/2013 based on a call with Joe.  The version
in this folder implements these concepts.


model model1.nam 1
model model2.nam 2
tdis model.tdis
sms model.sms 1
xrs model.xrs 1 2
begin solution_group 1
  mxiter 1
  sms 1
    model 1
    model 2
end
_____________________
model model1.nam 1
model model2.nam 2
tdis model.tdis
xrs model.xrs 1 2
sms model1.sms 1
sms model2.sms 2
begin solution_group 1
  mxiter 10
  sms 1
    model 1
  sms 2
    model 2
end
_____________________
explicit
TIMESTEP

  MODEL CONVERGENCE PICARD

    SMS1 (SWR)
      subtimesteps
        PICARD/NEWTON (kiter loop)
          formulate
          newton-raphson
          linsolve

    SMS2 (MODFLOW)
      subtimesteps
        PICARD/NEWTON (kiter loop)
          formulate
          newton-raphson
          linsolve

    CROSS

_____________________
implicit
TIMESTEP

  MODEL CONVERGENCE PICARD

    SMS1 (SWR, MODFLOW)
        PICARD/NEWTON (only 1 iteration)
          formulate
          newton-raphson
          linsolve

    CROSS


______________________
for kper=1,nper
  for kstp=1,nstp

    for sg in solution_group
      for kiter in xrange(sg.mxiter)
        for sms in sg.solutions
          for smsiter=1,smsmaxiter
            for model in sms.models
              model.formulate
              model.nr
            sms.linsolve
        for c in sg.crosses
          c.fill
