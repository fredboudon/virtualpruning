def penumerate(values, percentdelta = 1):
    nbelem = float(len(values))
    prevpercent = 0
    for current, value in enumerate(values):
        curpercent = 100. * (current + 1) / nbelem
        if curpercent - prevpercent >= percentdelta:
            print ("Done %.2f%s ..." % (curpercent, '%')) 
            prevpercent = curpercent
        yield current, value
