###############################################################################
# NAME:         standard_curves.py
# VERSION:      2.0 beta 3
# RELEASED:     21 SEPTEMBER 2007
# MODIFIED:     21 SEPTEMBER 2007
# AUTHORS:      John B. Cole (john.cole@ars.usda.gov)
# DESCRIPTION:  Plots the standard curve for a given breed-trait-DIM combination.
#               This program is part of the bestpred package from AIPL.
###############################################################################
import glob, math, string, sys

#-- Keep this block as-is or you'll break
#-- your graphs.
import matplotlib.numerix.ma as M
import matplotlib
matplotlib.use('Agg')
import pylab
#--

# Change this to 1 if you want graphs plotted
# with kilograms instead of pounds.
kgtolb = 1

# These parameters affect the appearance of your
# graphs. If you want larger typefaces on the
# graphs increase the *.*size items.
params = {'backend': 'Agg',
        'axes.labelsize': 18,
        'text.fontsize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14
        }

number2name = {1:'Milk',
               2:'Fat',
               3:'Protein',
               4:'SCS'
            }

number2breed = {1:'Ayrshire',
                2:'Brown Swiss',
                3:'Guernsey',
                4:'Holstein',
                5:'Jersey',
                6:'Milking Shorthorn'
            }

number2parity = {1:'first',
                 2:'later'
            }

number2units = {1:'lbs',
                0:'kgs'
            }

# Breed, trait, parity
woods_means =   {
                    1:{
                        1:{
                            1:{'a':14.7197,'b':.1312,'c':.0012},
                            2:{'a':29.9819,'b':-.0037,'c':.0008},
                        },
                        2:{
                            1:{'a':.5968,'b':.1053,'c':.0009},
                            2:{'a':1.2941,'b':-.0306,'c':.0007},
                        },
                        3:{
                            1:{'a':.3973,'b':.1619,'c':.0011},
                            2:{'a':.8529,'b':.0186,'c':.0008},
                        },
                    },
                    2:{
                        1:{
                            1:{'a':17.3838,'b':.1127,'c':.0009},
                            2:{'a':28.2813,'b':.0541,'c':.0010},
                        },
                        2:{
                            1:{'a':.7576,'b':.0862,'c':.0006},
                            2:{'a':1.2817,'b':.0149,'c':.0007},
                        },
                        3:{
                            1:{'a':.4743,'b':.1506,'c':.0008},
                            2:{'a':.8090,'b':.0787,'c':.0009},
                        },
                    },
                    3:{
                        1:{
                            1:{'a':18.7181,'b':.0614,'c':.0007},
                            2:{'a':28.8644,'b':.0085,'c':.0011},
                        },
                        2:{
                            1:{'a':.5704,'b':.1250,'c':.0006},
                            2:{'a':1.1929,'b':.0121,'c':.0009},
                        },
                        3:{
                            1:{'a':.4747,'b':.1085,'c':.0006},
                            2:{'a':.8454,'b':.0233,'c':.0009},
                        },
                    },
                    4:{
                        1:{
                            1:{'a':17.8808,'b':.1616,'c':.00107},
                            2:{'a':34.5682,'b':.0625,'c':.00126},
                        },
                        2:{
                            1:{'a':.7632,'b':.1473,'c':.00103},
                            2:{'a':1.4489,'b':.0516,'c':.00104},
                        },
                        3:{
                            1:{'a':.5538,'b':.1493,'c':.0010},
                            2:{'a':1.0863,'b':.0497,'c':.00124},
                        },
                    },
                    5:{
                        1:{
                            1:{'a':14.9604,'b':.1163,'c':.0009},
                            2:{'a':26.0002,'b':.0291,'c':.0010},
                        },
                        2:{
                            1:{'a':.5871,'b':.1410,'c':.0008},
                            2:{'a':1.1356,'b':.0280,'c':.0008},
                        },
                        3:{
                            1:{'a':.4149,'b':.1610,'c':.0008},
                            2:{'a':.9085,'b':.0225,'c':.0007},
                        },
                    },
                    6:{
                        1:{
                            1:{'a':29.0200,'b':.0116,'c':-.0001},
                            2:{'a':45.6060,'b':-.0508,'c':.0010},
                        },
                        2:{
                            1:{'a':1.0724,'b':.1884,'c':-.0001},
                            2:{'a':2.5053,'b':-.1680,'c':.0003},
                        },
                        3:{
                            1:{'a':.9473,'b':.0092,'c':.0002},
                            2:{'a':1.3777,'b':-.0524,'c':.0010},
                        },
                    },
                }

# Means and SD for SCS are modelled using curve C$
# from Morant and Gnanasankthy (19889).
# Breed, parity
mandg_means =   {
                    1:{
                        1:{'a':1.9724,'b':-.00287,'c':-1.42E-6,'d':7.9803},
                        2:{'a':1.8346,'b':-.00855,'c':-0.00002,'d':15.7619},
                    },
                    2:{
                        1:{'a':1.522,'b':-.00443,'c':-4.9E-6,'d':19.7959},
                        2:{'a':2.6844,'b':-.00481,'c':-6.58E-6,'d':15.5049}
                    },
                    3:{
                        1:{'a':2.3242,'b':-.00095,'c':6.322E-7,'d':16.6699},
                        2:{'a':2.8001,'b':-.00301,'c':-2.95E-6,'d':7.4281}
                    },
                    4:{
                        1:{'a':1.9178,'b':-.00253,'c':-1.36E-6,'d':18.5222},
                        2:{'a':3.0034,'b':-.00333,'c':-3.75E-6,'d':9.9904}
                    },
                    5:{
                        1:{'a':2.4263,'b':-.00178,'c':-3.5E-7,'d':16.3506},
                        2:{'a':2.8626,'b':-.00371,'c':-4.86E-6,'d':14.6406}
                    },
                    6:{
                        1:{'a':0.3585,'b':-.01150,'c':-0.00002,'d':38.6850},
                        2:{'a':0.6391,'b':-.00531,'c':-6.01E-6,'d':46.2029}
                    },
                }

# Simple command line processing to get files.
if len(sys.argv) < 5:
    print '[ERROR]: You must provide a breed (1-6), trait (1-4), parity (1,2), and DIM.'
    sys.exit(0)
else:
    breed = int(sys.argv[1])
    trait = int(sys.argv[2])
    parity = int(sys.argv[3])
    dimin = int(sys.argv[4])

if breed not in xrange(1,7):
    breed = 4 				    # Holstein
if trait not in xrange(1,5):
    trait = 1				    # MFPS
if parity not in xrange(1,3):
    parity = 1                   # Mean
if dimin < 0 or dimin > 999:
    dimin = 999

std, dim = [], []
plotstd, plotdim = [], []
stdsum = 0.

def compute_curve_point(breed,trait,parity,day):
    y = 0.
    if trait < 4:
        # Woods curve for MFP
        a = woods_means[breed][trait][parity]['a']
        b = woods_means[breed][trait][parity]['b']
        c = woods_means[breed][trait][parity]['c']
        y = a * day**b * math.exp(-day*c)
    else:
        # Morant and Gnanasakthy curve for SCS
        day10 = day + 10
        a = mandg_means[breed][parity]['a']
        b = mandg_means[breed][parity]['b']
        c = mandg_means[breed][parity]['c']
        d = mandg_means[breed][parity]['d']
        y = a - ( b*day10 ) + (( c*day10**2 ) / 2) + ( d/day10 )
    return y

if __name__ == '__main__':

    for day in xrange(1,dimin+1):
        dim.append(day)
        point = compute_curve_point(breed,trait,parity,day)
        std.append(point)
        stdsum = stdsum + point

    plotdim = M.array(dim)
    plotstd = M.array(std)
    maxdim = plotdim.max()

    if kgtolb == 1 and trait < 4:
        plotstd = plotstd * 2.2
        stdsum = stdsum * 2.2

    # Plot standard curves
    pylab.rcParams.update(params)
    fig = pylab.figure()
    plot = fig.add_subplot(111)
    plot_title = 'Standard curve for %s parity %s %s (%s d)' % ( number2parity[parity], \
        number2breed[breed], number2name[trait], dimin )
    pylab.title(plot_title)
    pylab.xlabel('DIM')
    if trait < 4:
        pylab.ylabel('%s Yield (%s)' % ( number2name[trait], number2units[kgtolb] ) )
    else:
        pylab.ylabel('%s' % ( number2name[trait] ) )
    plot.plot(plotdim, plotstd, '-g', linewidth=2)
    figfile = 'curve_std_%s_%s_%s_%s.png' % ( number2breed[breed], number2name[trait], \
        number2parity[parity], dimin )
    if trait < 4:
        # Automatically scale axes.
        pylab.axis([1,maxdim,0.,plotstd.max()+5])
        plot.annotate(str(round(stdsum)), xy=(0.85, 0.95),  xycoords='axes fraction', color='green', size=16)
    else:
        pylab.axis([1,maxdim,0.,plotstd.max()+0.1])
        plot.annotate(str(round(stdsum/dimin,2)), xy=(0.85, 0.85),  xycoords='axes fraction', color='green', size=16)
    print '\rWriting file %s' % ( figfile )
    pylab.savefig(figfile)
    pylab.show()