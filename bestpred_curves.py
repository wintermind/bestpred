###############################################################################
# NAME:         aiplcurves.py
# VERSION:      2.0 beta 3
# RELEASED:     01 AUGUST 2007
# MODIFIED:     10 September 2007
# AUTHORS:      John B. Cole (john.cole@ars.usda.gov)
# DESCRIPTION:  Reads output files from bestpred and plots the resulting
#               lactation curves. This program is part of the bestpred package
#               from AIPL.
###############################################################################
import glob, string, sys

#-- Keep this block as-is or you'll break
#-- your graphs.
import matplotlib.numerix.ma as M
import matplotlib
matplotlib.use('Agg')
import pylab
#--

# Change this to 1 if you want graphs plotted
# with kilograms instead of pounds.
lbtokg = 0

simmast = 1

# These parameters affect the appearance of your
# graphs. If you want larger typefaces on the
# graphs increase the *.*size items.
params = {'backend': 'Agg',
        'axes.labelsize': 18,
        'text.fontsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14
        }

number2name = {'1':'Milk',
                '2':'Fat',
                '3':'Protein',
                '4':'SCS'
            }

number2units = {1:'kgs',
                0:'lbs'
            }

#
# Simple command line processing to get files.
# !!! If you change the CURVEfile parameter in
#     bestpred.par you'll also need to change
#     "cowcurve" below if you don't pass in actual
# .   filenames and let aiplcurves glob the input.
#
if len(sys.argv) == 1:
    #print '[ERROR]: You must provide a plot style: [C]omplete or [P]artial'
    print '[ERROR]: You must provide a trait: 1 = milk, 2 = fat, 3 = protein, 4 = SCS, 5 = all'
    sys.exit(0)
elif len(sys.argv) == 2:
    curvefiles = glob.glob("cowcurve.*")
else:
    curvefiles = sys.argv[2:]

# matplotlib's line plots, created using the plot command,
# do not know what to do with masked arrays. So the plotstyle
# is used to select either line or scatter plots. This may be
# fixed in newer versions of matplotlib.
#if sys.argv[1] not in ['C','P']:
#    plotstyle = 'P'
#else:
#    plotstyle = sys.argv[1]
if sys.argv[1] not in ['1','2','3','4','5']:
    plottrait = 1
else:
    plottrait = sys.argv[1]

#
# Loop over the filenames provided and create a graph from each file.
#
threshold = -999.0
lastcowid = 0
lastlacnum = 0
lasttrait = 0
nround = 0
dim, dbp, tdy, std, lac, dev = [], [], [], [], [], []
stdsum, dbpsum = 0., 0.

for cf in curvefiles:
    print 'Processing %s' % ( cf )
    infile = open(cf,'r')
    while 1:
        line = infile.readline()
        if not line:
            break
        lpre = " ".join(line.split())
        lp = string.split(lpre)
        if len(lp) > 0:
            cowid = lp[0]
            lacnum = lp[1]
            trait = lp[2]
            if nround == 0:
                lasttrait = trait
                lastcowid = cowid
                lastlacnum = lacnum
                nround = nround + 1
            if cowid != lastcowid or trait != lasttrait:
                if int(lasttrait) == int(plottrait) or int(plottrait) == 5:
                    if cowid != lastcowid:
                        #print 'Cowid changed from %s to %s' % ( lastcowid, cowid )
                        _cowid = lastcowid
                        _lacnum = lastlacnum
                    else:
                        _cowid = cowid
                        _lacnum = lacnum
                    if trait != lasttrait:
                        #print 'Trait changed from %s to %s' % ( lasttrait, trait )
                        _trait = lasttrait
                    else:
                        _trait = trait
                    # We're processing a new animal
                    plotdim = M.array(dim)
                    plottd_ = M.array(tdy)
                    plottdy = M.masked_where(plottd_ == threshold, plottd_)
                    plotla_ = M.array(lac)
                    plotlac = M.masked_where(plotla_ == -2.0, plotla_)
                    plotst_ = M.array(std)
                    plotstd = M.masked_where(plotst_ == threshold, plotst_)
                    plotdb_ = M.array(dbp)
                    plotdbp = M.masked_where(plotdb_ == threshold, plotdb_)
                    plotde_ = M.array(dev)
                    plotdev = M.masked_where(plotde_ == threshold, plotde_)
                    maxdim = plotdim.max()

                    if lbtokg == 1 and int(_trait) < 4:
                        plottdy = plottdy / 2.2
                        plotlac = plotlac / 2.2
                        plotstd = plotstd / 2.2
                        stdsum = stdsum / 2.2
                        dbpsum = dbpsum / 2.2

                    # Plot lactation curves for each cow
                    # This will only work for lactation curve files with the
                    # form, e.g., 'cowcurve.S.HOUSA.EX.COW.0041.MT'.
                    #fig = pylab.figure(frameon=False)
                    pylab.rcParams.update(params)
                    fig = pylab.figure()
                    plot = fig.add_subplot(111)
                    plot_title = 'BP of %s for Cow %s (lactation %s)' % ( number2name[_trait], _cowid, _lacnum )
                    pylab.title(plot_title)
                    pylab.xlabel('DIM')
                    if int(_trait) < 4:
                        pylab.ylabel('%s Yield (%s)' % ( number2name[_trait], number2units[lbtokg] ) )
                    else:
                        pylab.ylabel('%s Yield' % ( number2name[_trait] ) )
                    plot.plot(plotdim, plottdy, 'go')
                    #if plotstyle == 'C':
                    plot.plot(plotdim, plotlac, '-g', linewidth=2)
                    plot.plot(plotdim, plotstd, '-b', linewidth=2)
                    ##plot.plot(plotdim, plotdev, '-r', linewidth=1)
                    ##plot.plot(plotdim, plotdbp, '-y', linewidth=1)
                    #else:
                    #    plot.plot(plotdim, plotlac, 'o-g')
                    #    plot.plot(plotdim, plotstd, 'o-b')
                    # Automatically scale axes.
                    pylab.axis([1,maxdim,0.,max(plotlac.max(),plotstd.max(),plottdy.max())])
                    figfile = 'curve_%s_%s_%s.png' % ( _cowid,number2name[_trait],_lacnum )
                    if int(_trait) < 4:
                        plot.annotate(str(round(dbpsum)), xy=(0.85, 0.95),  xycoords='axes fraction', color='green', size=16)
                        plot.annotate(str(round(stdsum)), xy=(0.85, 0.90),  xycoords='axes fraction', color='blue', size=16)
                    else:
                        plot.annotate(str(round(dbpsum/365.,2)), xy=(0.85, 0.85),  xycoords='axes fraction', color='green', size=16)
                        plot.annotate(str(round(stdsum/365.,2)), xy=(0.85, 0.80),  xycoords='axes fraction', color='blue', size=16)
                    print '\rWriting file %s' % ( figfile )
                    pylab.savefig(figfile)
                    pylab.show()

                # We're not going to plot this trait
                else:
                    pass

                if trait != lasttrait:
                    lasttrait = trait
                if cowid != lastcowid:
                    lastcowid = cowid
                    lastlacnum = lacnum

                # We still need to clean up
                if nround > 1:
                    # Cleanup to avoid multiple-plots-on-the-same-canvas
                    # problem.
                    try: del plotdim
                    except: pass
                    try: del plottd_
                    except: pass
                    try: del plottdy
                    except: pass
                    try: del plotlac
                    except: pass
                    try: del plotstd
                    except: pass
                    try: del plotdbp
                    except: pass
                    try: del plotde_
                    except: pass
                    try: del plotdev
                    except: pass
                    try: del plotdb_
                    except: pass
                    try: del plotla_
                    except: pass
                    try: del plotst_
                    except: pass
                    try: del dim
                    except: pass
                    try: del tdy
                    except: pass
                    try: del lac
                    except: pass
                    try: del std
                    except: pass
                    try: del dbp
                    except: pass
                    try: del dev
                    except: pass
                    try: del dev
                    except: pass
                    try: del fig
                    except: pass
                    try: del stdsum
                    except: pass
                    try: del dbpsum
                    except: pass

                    dim, dbp, tdy, std, lac, dev = [], [], [], [], [], []
                    stdsum, dbpsum = 0., 0.

            else:
                # We're adding another line of data to an existing
                # cow-trait combination.
                dim.append(int(lp[3]))
                dbp.append(float(lp[4]))
                tdy.append(float(lp[5]))
                std.append(float(lp[6]))
                lac.append(float(lp[4])+float(lp[6]))
                dev.append(float(lp[7]))
                # Increment sums
                if int(lp[3]) <= 305:
                    stdsum = stdsum + float(lp[6])
                    dbpsum = dbpsum + ( float(lp[4]) + float(lp[6]) )
            nround = nround + 1
        else:
            pass
    infile.close()

    # We ned to catch the data from the last animal in the input file
    # and plot it. Wow, this is ugly.
    plotdim = M.array(dim)
    plottd_ = M.array(tdy)
    plottdy = M.masked_where(plottd_ == threshold, plottd_)
    plotla_ = M.array(lac)
    plotlac = M.masked_where(plotla_ == -2.0, plotla_)
    plotst_ = M.array(std)
    plotstd = M.masked_where(plotst_ == threshold, plotst_)
    plotdb_ = M.array(dbp)
    plotdbp = M.masked_where(plotdb_ == threshold, plotdb_)
    plotde_ = M.array(dev)
    plotdev = M.masked_where(plotde_ == threshold, plotde_)
    maxdim = plotdim.max()

    if lbtokg == 1 and int(_trait) < 4:
        plottdy = plottdy / 2.2
        plotlac = plotlac / 2.2
        plotstd = plotstd / 2.2
        stdsum = stdsum / 2.2
        dbpsum = dbpsum / 2.2

    pylab.rcParams.update(params)
    fig = pylab.figure()
    plot = fig.add_subplot(111)
    plot_title = 'BP of %s for Cow %s (lactation %s)' % ( number2name[trait], cowid, lacnum )
    pylab.title(plot_title)
    pylab.xlabel('DIM')
    if int(trait) < 4:
        pylab.ylabel('%s Yield (%s)' % ( number2name[trait], number2units[lbtokg] ) )
    else:
        pylab.ylabel('%s' % ( number2name[trait] ) )
    plot.plot(plotdim, plottdy, 'go')
    #if plotstyle == 'C':
    plot.plot(plotdim, plotlac, '-g', linewidth=2)
    plot.plot(plotdim, plotstd, '-b', linewidth=2)
    #plot.plot(plotdim, plotdev, '-r', linewidth=1)
    #plot.plot(plotdim, plotdbp, '-y', linewidth=1)
    #else:
    #    plot.plot(plotdim, plotlac, 'o-g')
    #    plot.plot(plotdim, plotstd, 'o-b')
    # Automatically scale axes.
    pylab.axis([1,maxdim,0.,max(plotlac.max(),plotstd.max(),plottdy.max())])
    figfile = 'curve_%s_%s_%s.png' % ( cowid,number2name[trait],lacnum )
    if int(_trait) == 4:
        plot.annotate(str(round(dbpsum/365.,2)), xy=(0.85, 0.85),  xycoords='axes fraction', color='green', size=16)
        plot.annotate(str(round(stdsum/365.,2)), xy=(0.85, 0.80),  xycoords='axes fraction', color='blue', size=16)
    else:
        plot.annotate(str(round(dbpsum)), xy=(0.85, 0.95),  xycoords='axes fraction', color='green', size=16)
        plot.annotate(str(round(stdsum)), xy=(0.85, 0.90),  xycoords='axes fraction', color='blue', size=16)
    print '\rWriting file %s' % ( figfile )
    pylab.savefig(figfile)
    pylab.show()
