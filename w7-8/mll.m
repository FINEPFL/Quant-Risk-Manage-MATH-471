function ret = mll(pse, x, copulaname)
    ret = -sum(log(copulapdf(copulaname, pse, x)));
end