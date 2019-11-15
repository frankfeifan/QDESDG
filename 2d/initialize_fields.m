function ssd = initialize_fields(w, time, vel, tau, iter, time_step, pre)

if nargin < 7
  pre = '';
end

%%% FAULT VARIABLES %%%
fn = [pre,'w.dat'];
ssd.w = SaveStreamData('Init',fn);
SaveStreamData('Write',ssd.w,w);
gn = [pre,'Time.dat'];
ssd.time = SaveStreamData('Init',gn);
SaveStreamData('Write',ssd.time,time);
hn = [pre,'Vel.dat'];
ssd.vel = SaveStreamData('Init',hn);
SaveStreamData('Write',ssd.vel,vel);

jn = [pre,'tau.dat'];
ssd.tau = SaveStreamData('Init',jn);
SaveStreamData('Write',ssd.tau,tau);

vvn = [pre,'iterations.dat'];
ssd.iterations= SaveStreamData('Init',vvn);
SaveStreamData('Write',ssd.iterations,iter);

dtn = [pre,'time_step.dat'];
ssd.time_step= SaveStreamData('Init',dtn);
SaveStreamData('Write',ssd.time_step,time_step);

