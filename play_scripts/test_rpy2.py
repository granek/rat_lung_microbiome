import rpy2.robjects as robjects

print(robjects.r['R.home']()[0])
