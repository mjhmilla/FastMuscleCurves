function  dsee = calcDseeUmat41(q,cm,A_rel,B_rel,F_SEE,F_PEE,F_isom,
     1 dot_l_MTC,D0,C2,C1,C0,d_SE_max)


%        real*8 q,A_rel,B_rel,F_SEE,F_PEE,F_isom,
%     1   dot_l_MTC,d_SE_max,D0,C2,C1,C0
%    Force-dependent serial (SE)-Damping

D0 = cm(10)*B_rel*d_SE_max*(cm(28)+(1.0-cm(28))*(q*F_isom+F_PEE/cm(9)))
C2 = d_SE_max*(cm(28)-(A_rel-F_PEE/cm(9))*(1.0-cm(28)))
C1 = -C2*dot_l_MTC-D0-F_SEE+F_PEE-cm(9)*A_rel
C0 = D0*dot_l_MTC+cm(10)*B_rel*(F_SEE-F_PEE-cm(9)*q*F_isom)
