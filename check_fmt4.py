# Read a "reduced" format 4 record, split it up, and write it to STDOUT
# for checking BESTPRED inputs.

if __name__ == '__main__':

	infile = 'format4.dat'
	ifh = open(infile, 'r')
	linecount = 0
	for line in ifh:
		linecount += 1
		print 'Line %s is of length %s' % ( linecount, len(line) )
		print '\t17-byte ID  : "%s"' % ( line[2:19] )
		print '\tBirth date  : "%s"' % ( line[70:78] )
		print '\tHerd code   : "%s"' % ( line[106:114] )
		print '\tCalving date: "%s"' % ( line[127:135] )
		print '\tDIM         : "%s"' % ( line[135:138] )
		print '\tLactation   : "%s"' % ( line[158:160] )
		print '\tPrevious DO : "%s"' % ( line[245:248] )
		print '\tNumber of TD: "%s"' % ( line[248:250] )
		numtd = int(line[248:250])
		for curtd in xrange(0,numtd):
			print '\tTest day segment %s' % ( curtd )
			offset = 250 + ( curtd * 23 )
			print '\t\tSegment           : "%s"' % ( line[offset:offset+23] )
			print '\t\tTD DIM            : "%s"' % ( line[offset:offset+3] )
			print '\t\tSupervision code  : "%s"' % ( line[offset+3:offset+4] )
			print '\t\tMilking frequency : "%s"' % ( line[offset+5:offset+6] )
			print '\t\tMilkings weighed  : "%s"' % ( line[offset+6:offset+7] )
			print '\t\tMilking sampled   : "%s"' % ( line[offset+7:offset+8] )
			print '\t\tMilk-recorded days: "%s"' % ( line[offset+8:offset+10] )
			print '\t\tActual milk       : "%s"' % ( line[offset+13:offset+17] )
			print '\t\tActual fat pct    : "%s"' % ( line[offset+17:offset+19] )
			print '\t\tActual protein pct: "%s"' % ( line[offset+19:offset+21] )
			print '\t\tActual SCS        : "%s"' % ( line[offset+21:offset+23] )
