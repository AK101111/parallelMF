def timetosec(t):
    t = t.split('m')
    return (float(t[0])*60) + float(t[1][:-2])

with open('pp.txt','r') as file:
    x = 1
    d = []
    l = []
    for line in file:
    	if x == 1:
    	    #print(line.split())
    	    error = float(line.split()[2])
    	    iter = int(line.split()[4])
    	    x += 1
    	elif x == 2:
    		#print(line.split())
    		x += 1
    		continue
    	elif x == 3:
    		#print(line.split())
    		t = line.split()[1]
    		x += 1
    	elif x == 4:
    		#print(line.split())
    		x += 1
    		continue
    	elif x == 5:
    		#print(line.split())
    		d = [t,iter, error]
    		x = 1
    		l.append(d)
    for x in l:
    	print ' '.join(str(xx)+',' for xx in x)