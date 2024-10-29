function dataHL2002 = getHerzogLeonard2002Keypoints()

lceOpt = 42.8571; %The optimal fiber length is indirectly reported in
                  % Herzog and Leonard 2002
                  %27mm/s is listed as 63\% of 1 optimal fiber length/s
                  %in Fig 3 caption, and 9mm is listed as 21% of optimal
                  %fiber length: 9/0.21 = 27/0.63 = 42.8571 mm

lstartN = 0.9607448010823296; %Found by hand: this fits the force-length curves
                              %much better. It is plausible that there is
                              % 5% error in the quoted value of lceOpt (2mm)
                              % both because only a few points were taken, and
                              % also because the exact value depends on how
                              % this point is defined. Beause the force length
                              % relation has a plateau, there is some abiguity
                              % of where to put the optimal fiber length.
                              % I tend to put it in the middle of the plateau
                              % while others put it at the start.

fceOpt = 41.661837;

%From the text
dataHL2002.l = [0,3,6,9];


dataHL2002.lce    = ([0,3,6,9]+lstartN*lceOpt)./lceOpt;

dataHL2002.dl    = diff(dataHL2002.lce);

%Manually digitized
dataHL2002.fpe  = [mean(0.422,0.3177,0.3769),...
                   mean(1.0499,1.0377,1.0934),...
                   mean(3.0424,2.774,2.8157),...
                   mean(7.6554,7.2932,7.4896)];

dataHL2002.fmt  = [mean(21.9865,21.4345,21.0701),...
                   mean(27.0183,27.0778,26.2587),...
                   mean(25.7527,25.4432,24.9761),...
                   23.9485];



dataHL2002.fa   = dataHL2002.fmt-dataHL2002.fpe;

dataHL2002.fpeN =dataHL2002.fpe./fceOpt;
dataHL2002.fmtN =dataHL2002.fmt./fceOpt;
dataHL2002.faN  =dataHL2002.fa ./fceOpt;

