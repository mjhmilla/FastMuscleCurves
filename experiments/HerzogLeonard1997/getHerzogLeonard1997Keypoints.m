function dataHL1997 = getHerzogLeonard1997Keypoints()

dataHL1997.l    = [0.688979,0.9011];
dataHL1997.dl   = dataHL1997.l-0.688979;
dataHL1997.fpe  = [0.0237796,0.0760336];
dataHL1997.fa   = [0.770336,1.03486]-dataHL1997.fpe;
dataHL1997.dfa   = dataHL1997.fa-dataHL1997.fa(1,1);
dataHL1997.rectangle.fl =   [dataHL1997.dl(1,1),dataHL1997.fa(1,1);...
                             dataHL1997.dl(1,2),dataHL1997.fa(1,1);...
                             dataHL1997.dl(1,2),dataHL1997.fa(1,2);...
                             dataHL1997.dl(1,1),dataHL1997.fa(1,2);...
                             dataHL1997.dl(1,1),dataHL1997.fa(1,1)];
dataHL1997.rectangle.fpe =  [dataHL1997.dl(1,1),dataHL1997.fpe(1,1);...
                             dataHL1997.dl(1,2),dataHL1997.fpe(1,1);...
                             dataHL1997.dl(1,2),dataHL1997.fpe(1,2);...
                             dataHL1997.dl(1,1),dataHL1997.fpe(1,2);...
                             dataHL1997.dl(1,1),dataHL1997.fpe(1,1)];