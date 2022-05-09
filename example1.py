import coaxcal.coax_filter as cf

myCf = cf.CoaxFilter(11.5E-3, 1E-3,12E-3, 2.25, 1)

print("### Dataset ###")
print("Internal cavity: ", myCf.cavRad * 1E3, " mm")
print("Capacitance per meter (at maximum radius ", myCf.maxRad * 1E3, " mm): ", myCf.getCap(1) * 1E9, " nF")
print("Inductance per meter (at minimum radius ", myCf.minRad * 1E3, " mm): ", myCf.getInd(1) * 1E9, " nH")

myFilter = [('L' , 15.21e-9), ('C', 7.275E-12), ('L' , 26.19e-9), ('C', 7.275E-12), ('L' , 15.21e-9)]
myCf.setLPF(myFilter)

myCoaxFilter = []

print("### LC Filter structure ###")
for i in myCf.lpfList:
    if i[0] == 'L':
        curInd = myCf.getIndLen(i[1])
        print("Inductance value: ", i[1] * 1e9, " nH, Length: ", curInd * 1E3, " mm")
        myCoaxFilter.append(('L', curInd))
    elif i[0] == 'C':
        curCap = myCf.getCapLen(i[1])
        print("Capacitance value: ", i[1] * 1e9, " nF, Length: ", curCap * 1E3, " mm")
        myCoaxFilter.append(('L', curCap))

myCf.setCoaxLPF(myCoaxFilter)