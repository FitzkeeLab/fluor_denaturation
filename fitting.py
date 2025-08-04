#
# title:   fitting.py
# summary: Model Fitting Script
# author:  Nick Fitzkee (nfitzkee at chemistry.msstate.edu)
# date:    February 17, 2021
#

import numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
import scipy.optimize as opt
import random, math

R = 1.9872e-3    # Gas constant in kcal mol-1 K-1
T = 273.15+22.0  # Temperature in K
RT = R*T

def load_data(infile, fixed_err=None, frac_err=None, guess=0.05):

    close_f = False

    if type(infile) is type(''):
        f = open(infile)
        close_f = True
    else:
        f = infile

    l = f.readline()

    result = []
    avg = 0.0
    ymin = None
    col3 = 0
    col2 = 0

    while l:
        l = l.strip()

        if not l or l[0] == '#':
            l = f.readline()
            continue

        toks = l.split('#')[0]
        toks = toks.split()

        if len(toks) == 2:
            col2 = col2 + 1
            x, y = float(toks[0]), float(toks[1])
            err = 0.0

            if not frac_err is None:
                err = frac_err*y
            elif not fixed_err is None:
                err = fixed_err

        elif len(toks) == 3:
            col3 = col3 + 1
            x, y, err = float(toks[0]), float(toks[1]), float(toks[2])
        else:
            raise BaseException("Unexpected column count!")

        if ymin is None or abs(y) < ymin:
            ymin = abs(y)

        result.append((x, y, err))

        l = f.readline()

    if close_f:
        f.close()

    if col2 > 0 and col3 > 0:
        raise BaseException("Multiple column counts detected!")

    if col2 > 0 and frac_err is None and fixed_err is None:
        guess_err = ymin * guess

        for i in range(len(result)):
            x, y, err = result[i]
            result[i] = (x, y, guess_err)

    return result

def boot_resample_data(orig_pairs):
    result = []

    for i in range(len(orig_pairs)):
        pair = random.choice(orig_pairs)
        result.append(pair)

    return result

def err_resample_data(orig_pairs):
    result = []

    for i in range(len(orig_pairs)):
        x, y, err = orig_pairs[i]
        y_new = random.normalvariate(y, err)
        result.append((x, y_new, err))

    return result

def fy(x, par, idx=0):
    af, bf, au, bu, dGH2O, m = par

    dG = dGH2O - m * x
    Keq = math.exp(-dG/RT)
    ff = 1.0/(1.0+Keq)
    fu = Keq/(1.0+Keq)

    return (af*x+bf)*ff + (au*x+bu)*fu

def chi_sqr(par, data):
    return sum([((y-fy(x, par))/w)**2 for (x, y, w) in data])

def plot_residuals(data, par, inter=False, output=None, cfname=None):
    
    plt.close()

    resid = [y-fy(x, par) for (x, y, w) in data]
    x, y, err = zip(*data)

    if cfname:
        with open(cfname, 'w') as cfile:
            for i in range(len(x)):
                cfile.write('%13.5e %13.5e %13.5e\n' %
                            (x[i], resid[i], err[i]))
    
    fig = plt.figure()
    plt.errorbar(x, resid, yerr=err, fmt='.') 

    if inter:
        plt.show()

    if output:
        plt.savefig(output)

    plt.close(fig)
    
def plot_transform(data, par, inter=None, plot=None, points=None, fit=None):

    plt.close()
    
    resid = [y-fy(x, par) for (x, y, w) in data]
    ptx, orgpty, orgerr = zip(*data)

    new_par = list(par)[:]

    new_par[0] = 0.0
    new_par[1] = 1.0
    new_par[2] = 0.0
    new_par[3] = 0.0
    
    trnpty = []
    trnerr = []

    if points:
        pfile = open(points, 'w')

    ymin = min(orgpty)
    ymax = max(orgpty)
    scale = 1.0 / (ymax - ymin)
    
    for i in range(len(ptx)):
        #print('original: %10.3f %10.3f %10.3f %10.3f' %
        #      (ptx[i], orgpty[i], orgerr[i], resid[i]))
        
        new_y = fy(ptx[i], new_par)+resid[i]*scale
        new_err = orgerr[i]*scale

        trnpty.append(new_y)
        trnerr.append(new_err)

        #print('transfrm: %10.3f %10.3f %10.3f %10.3f' %
        #      (ptx[i], new_y, new_err, resid[i]*scale)) 

        if points:
            pfile.write('%13.5e %13.5e %13.5e\n' % (ptx[i], new_y, new_err))

    if points:
        pfile.close()
        
    if fit:
        ffile = open(fit, 'w')
        
    curve = []
    xmin  = 0.0
    xmax  = max(ptx)*1.05
    dx    = 0.01
    x     = xmin

    while x <= xmax:
        y = fy(x, new_par)
        curve.append((x, y))

        if fit:
            ffile.write('%13.5e %13.5e\n' % (x, y))
            
        x = x+dx

    if fit:
        ffile.close()

    lnx, lny = zip(*curve)
    
    fig = plt.figure()
    plt.errorbar(ptx, trnpty, yerr=trnerr, fmt='.') 
    plt.plot(lnx, lny)

    if inter:
        plt.show()
        
    if plot:
        plt.savefig(plot)

    plt.close(fig)



def plot_fit(data, par, inter=False, output=None, cfname=None):
    ptx, pty, pterr = zip(*data)

    plt.close()
    
    cfile = None
    curve = []
    xmin  = 0.0
    xmax  = max(ptx)*1.05
    dx    = 0.01
    x     = xmin

    if cfname:
        cfile = open(cfname, 'w')
    
    while x <= xmax:
        y = fy(x, par)
        curve.append((x, y))

        if cfname:
            cfile.write('%13.5e %13.5e\n' % (x, y))
            
        x = x+dx

    if cfname:
        cfile.close()
        
    lnx, lny = zip(*curve)
    
    fig = plt.figure()
    plt.errorbar(ptx, pty, yerr=pterr, fmt='.') 
    plt.plot(lnx, lny)

    if inter:
        plt.show()

    if output:
        plt.savefig(output)

    plt.close(fig)    


def collect_parameter_stats(plist, pbest, pnames, inter=False, output=False,
                            fnames=False):
    pcols = list(zip(*plist))

    assert(len(pcols) == len(pnames))

    with open('par_summary.txt', 'w') as psf:
        psf.write('# %3s %-13s %-13s %-13s %-9s %-13s %-13s\n' %
                  ('var', 'best', 'avg', 'stddev', 'pcterr', '5%', '95%'))
    print('\n  %3s %-13s %-13s %-13s %-9s %-13s %-13s' %
          ('var', 'best', 'avg', 'stddev', 'pcterr', '5%', '95%'))
   
    for pidx in range(len(pnames)):
        pname = pnames[pidx]

        if fnames:
            with open('par_%s_list.txt' % pname, 'w') as f:
                for par in pcols[pidx]:
                    f.write('%13.5e\n' % par)

        vals = list(pcols[pidx][:])
        vals.sort()

        N = len(vals)
        avg = sum(vals)/N
        sd  = math.sqrt(sum([(x - avg)**2 for x in vals])/(N-1.0))
        ci5 = vals[int(0.05*N)]
        ci95 = vals[int(0.95*N)]

        with open('par_summary.txt', 'a') as psf:
            psf.write('%5s %13f %13f %13f %8.2f%% %13f %13f\n' %
                      (pname, pbest[pidx], avg, sd,
                       abs(sd/pbest[pidx]*100.0),
                       ci5, ci95))

        print('%5s %13f %13f %13f %8.2f%% %13f %13f' %
              (pname, pbest[pidx], avg, sd, abs(sd/pbest[pidx]*100.0),
               ci5, ci95))
        
        fig = plt.figure()
        plt.hist(pcols[pidx], 25)

        if output:
            plt.savefig('par_%s_histogram.pdf' % pname)
    
        if inter:
            plt.show()

def reset_error(data, new):
    result = []
    for x, y, err in data:
        result.append( (x, y, new) )

    return result

def main():
    data_file = "R2ab-Y722A-Den.txt"
    n_bootstrap = 100
    verbose = False

    pnames     = ['af',      'bf',  'au',    'bu', 'dG',  'm']
    parameters = [-0.25, 1.02, -0.001, 0.3, 3.0,  4.0]
    
    data = load_data(data_file)
    data = reset_error(data, 0.008)
    ndof = len(data) - len(parameters)
    
    print('* Starting Initial Fit...')

    plot_fit(data, parameters, True, 'fit.pdf', 'fit_curve.txt')
    print('  - Initial chi-sqr: %.3f' % chi_sqr(parameters, data))
    optres = opt.minimize(chi_sqr, parameters, args=data,
                        method='Nelder-Mead')

    print('  - Optimization Output:')
    print(optres)
    print('')
    
    print('  - Initial chi-sqr: %.3f' % chi_sqr(parameters, data))
    print('  - Final   chi-sqr: %.3f' % chi_sqr(optres.x, data))
    print('  - Deg. of Freedom: %i'   % ndof)
    print('  - Reduced chi-sqr: %.3f' % (chi_sqr(optres.x, data)/ndof))
    
    plot_fit(data, optres.x, False, 'fit.pdf', 'fit_curve.txt')
    plot_residuals(data, optres.x, False, 'residuals.pdf', 'residuals.txt')
    plot_transform(data, optres.x, inter=False, plot='transform.pdf',
                   points='trans_data.txt', fit='trans_curve.txt')

    ###########################
    # Begin Bootstrap Code

    par_list = []
    nb = 0
    nfail = 0
    
    print ('\n* Starting bootstrap (%i samples)...' % n_bootstrap)

    while nb < n_bootstrap:
        temp_data = boot_resample_data(data)

        b_optres = opt.minimize(chi_sqr, optres.x, args=temp_data,
                                method='Nelder-Mead')
        if b_optres.success:
            par_list.append(b_optres.x)
            nb = nb + 1
        else:
            if verbose:
                print('*** FAILED RESULT ***')
                print(b_optres)
            nfail = nfail + 1
    
    print('  - Failed fits: %i' % nfail)

    print('\n* Summary of parameters:')
    
    collect_parameter_stats(par_list, optres.x, pnames, inter=False,
                            output=True, fnames=True)

    
if __name__ == '__main__':
    main()
