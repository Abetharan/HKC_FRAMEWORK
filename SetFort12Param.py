def createFort12String(times):
    timeArray = [len(times)] + times
    strTimeArray = [str(i) for i in timeArray]
    strTimes = """ {} """.format("\n".join(strTimeArray))
    return(strTimes)
    