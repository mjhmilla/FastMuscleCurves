function y = calcDeGrooteFregly2016TendonMultiplier(ltN,m_kT,c1,c2,c3, der)

y = nan;

switch der
    case -1
        y = ((c1 * exp(m_kT * (ltN - c2))) / m_kT) - c3 * ltN;
        c = ((c1 * exp(m_kT * (1 - c2))) / m_kT) - c3 * 1;
        y = y - c;
    case 0
        y = c1 * exp(m_kT * (ltN - c2)) - c3;
    case 1
        y = c1 * m_kT * exp(m_kT * (ltN - c2));
    otherwise
        assert(0, 'Error: der must be -1, 0, 1');
end
